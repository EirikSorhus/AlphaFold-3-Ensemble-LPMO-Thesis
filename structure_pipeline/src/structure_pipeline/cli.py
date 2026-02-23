"""Command-line interface for the structure prediction pipeline."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .cases import (
    Case,
    CaseStatus,
    generate_cases,
    group_cases_by_ligand_model,
    group_cases_by_protein_model,
    load_cases_manifest,
    write_cases_manifest,
)
from .config import PipelineConfig, load_config
from .executors import LocalExecutor, SlurmExecutor
from .manifest.ligands import (
    LigandRecord,
    load_ccd_list,
    load_ligands_manifest,
    scan_ligands,
    write_ligands_manifest,
)
from .manifest.msa import (
    MSARecord,
    load_msa_manifest,
    match_msa,
    write_msa_manifest,
)
from .manifest.proteins import (
    ProteinRecord,
    load_proteins_manifest,
    parse_fasta,
    write_proteins_manifest,
)
from .resume import StatusTracker
from .runners import AF3Runner, BoltzRunner, RF3Runner
from .runners.mounts import compute_bind_mounts, log_mount_plan
from .runners.oligo import OligoRegistry

app = typer.Typer(
    name="structure-pipeline",
    help="Multi-model structure prediction pipeline for AF3, Boltz-2, and RF3",
    add_completion=False,
)
console = Console()


def get_config(config_path: Path | None) -> PipelineConfig:
    """Load configuration with error handling."""
    try:
        return load_config(config_path)
    except FileNotFoundError as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Configuration error:[/red] {e}")
        raise typer.Exit(1)


@app.command()
def version():
    """Show version information."""
    console.print(f"structure-pipeline v{__version__}")


@app.command()
def manifest(
    config: Annotated[
        Optional[Path],
        typer.Option("--config", "-c", help="Path to pipeline config YAML"),
    ] = None,
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing manifests"),
    ] = False,
    ligand_ccd_list: Annotated[
        Optional[Path],
        typer.Option(
            "--ligand-ccd-list",
            help=(
                "Path to a .txt file with one CCD code per line. "
                "When provided, ligands are taken from this list instead "
                "of scanning CIF files. All models will use standard CCD "
                "codes (no custom CIF paths)."
            ),
        ),
    ] = None,
):
    """Generate manifest files from input data.

    Parses FASTA, scans ligands, matches MSAs, and generates cases.
    """
    cfg = get_config(config)
    manifest_dir = cfg.outputs.manifest_dir

    console.print("[bold]Generating manifests...[/bold]")

    # 1. Parse proteins
    proteins_csv = manifest_dir / "proteins.csv"
    if proteins_csv.exists() and not force:
        console.print(f"  [yellow]Skipping[/yellow] proteins.csv (exists)")
        proteins = load_proteins_manifest(proteins_csv)
    else:
        console.print(f"  Parsing {cfg.inputs.proteins_fasta}...")
        proteins = list(parse_fasta(cfg.inputs.proteins_fasta))
        write_proteins_manifest(proteins, proteins_csv)
        console.print(f"  [green]✓[/green] {len(proteins)} proteins → {proteins_csv}")

    # 2. Scan ligands (or load from CCD list)
    ligands_csv = manifest_dir / "ligands.csv"
    if ligands_csv.exists() and not force:
        console.print(f"  [yellow]Skipping[/yellow] ligands.csv (exists)")
        ligands = load_ligands_manifest(ligands_csv)
    elif ligand_ccd_list is not None:
        console.print(f"  Loading CCD codes from {ligand_ccd_list}...")
        ligands = load_ccd_list(ligand_ccd_list)
        write_ligands_manifest(ligands, ligands_csv)
        console.print(
            f"  [green]✓[/green] {len(ligands)} ligands (CCD-list mode) → {ligands_csv}"
        )
    else:
        console.print(f"  Scanning {cfg.inputs.ligands_dir}...")
        ligands = list(scan_ligands(cfg.inputs.ligands_dir))
        write_ligands_manifest(ligands, ligands_csv)
        console.print(f"  [green]✓[/green] {len(ligands)} ligands → {ligands_csv}")

    # Show ligand categories
    categories = {}
    for lig in ligands:
        categories[lig.category] = categories.get(lig.category, 0) + 1
    for cat, count in sorted(categories.items()):
        console.print(f"    - {cat}: {count} ligands")

    # 3. Match MSAs (from directory or squashfs)
    msa_csv = manifest_dir / "msa.csv"
    if msa_csv.exists() and not force:
        console.print(f"  [yellow]Skipping[/yellow] msa.csv (exists)")
        msa_records = load_msa_manifest(msa_csv)
        missing_msa = [
            p.protein_id
            for p in proteins
            if p.protein_id not in {m.protein_id for m in msa_records}
        ]
    else:
        # Determine MSA sources
        msa_dir = cfg.inputs.msa_dir
        msa_sqsh_files = cfg.inputs.msa_sqsh_files

        if msa_sqsh_files:
            console.print(f"  Matching MSAs from {len(msa_sqsh_files)} squashfs file(s)...")
            for sqsh in msa_sqsh_files:
                console.print(f"    - {sqsh.name}")
        elif msa_dir:
            console.print(f"  Matching MSAs from {msa_dir}...")
        else:
            console.print("  [yellow]No MSA source configured[/yellow]")
            msa_records = []
            missing_msa = [p.protein_id for p in proteins]

        if msa_sqsh_files or msa_dir:
            msa_records, missing_msa = match_msa(
                proteins,
                msa_dir=msa_dir,
                msa_sqsh_files=msa_sqsh_files,
            )
            write_msa_manifest(msa_records, msa_csv)
            console.print(f"  [green]✓[/green] {len(msa_records)} MSAs matched → {msa_csv}")

    if missing_msa:
        console.print(
            f"  [yellow]Warning:[/yellow] {len(missing_msa)} proteins without MSA "
            "(will skip boltz/rf3 for these)"
        )
        for pid in missing_msa[:5]:
            console.print(f"    - {pid}")
        if len(missing_msa) > 5:
            console.print(f"    ... and {len(missing_msa) - 5} more")

    # 4. Generate cases
    cases_csv = manifest_dir / "cases.csv"
    if cases_csv.exists() and not force:
        console.print(f"  [yellow]Skipping[/yellow] cases.csv (exists)")
        cases = load_cases_manifest(cases_csv)
    else:
        console.print(f"  Generating cases for models: {cfg.models_enabled}...")
        cases = list(
            generate_cases(
                proteins=proteins,
                ligands=ligands,
                msa_records=msa_records,
                models=cfg.models_enabled,
                missing_msa_proteins=missing_msa,
            )
        )
        write_cases_manifest(cases, cases_csv)
        console.print(f"  [green]✓[/green] {len(cases)} cases → {cases_csv}")

    # Summary table
    table = Table(title="Case Summary")
    table.add_column("Model")
    table.add_column("Pending", justify="right")
    table.add_column("Skipped", justify="right")

    for model in cfg.models_enabled:
        model_cases = [c for c in cases if c.model == model]
        pending = sum(1 for c in model_cases if c.status == CaseStatus.PENDING.value)
        skipped = sum(1 for c in model_cases if c.status == CaseStatus.SKIPPED.value)
        table.add_row(model, str(pending), str(skipped))

    console.print(table)


@app.command()
def validate(
    config: Annotated[
        Optional[Path],
        typer.Option("--config", "-c", help="Path to pipeline config YAML"),
    ] = None,
):
    """Validate configuration and check all required paths."""
    cfg = get_config(config)

    console.print("[bold]Validating configuration...[/bold]\n")

    all_valid = True

    # Check input paths
    console.print("[bold]Input paths:[/bold]")

    # Required paths
    for name, path in [
        ("proteins_fasta", cfg.inputs.proteins_fasta),
        ("ligands_dir", cfg.inputs.ligands_dir),
    ]:
        exists = path.exists()
        status = "[green]✓[/green]" if exists else "[red]✗[/red]"
        console.print(f"  {status} {name}: {path}")
        if not exists:
            all_valid = False

    # Optional: msa_dir (can be None if using squashfs)
    if cfg.inputs.msa_dir:
        exists = cfg.inputs.msa_dir.exists()
        status = "[green]✓[/green]" if exists else "[red]✗[/red]"
        console.print(f"  {status} msa_dir: {cfg.inputs.msa_dir}")
        if not exists:
            all_valid = False

    # Check squashfs files
    for sqsh in cfg.inputs.msa_sqsh_files:
        exists = sqsh.exists()
        status = "[green]✓[/green]" if exists else "[red]✗[/red]"
        console.print(f"  {status} msa_sqsh: {sqsh}")
        if not exists:
            all_valid = False

    # Check model paths
    console.print("\n[bold]Model paths:[/bold]")

    if "af3" in cfg.models_enabled:
        runner = AF3Runner(cfg)
        errors = runner.validate_paths()
        if errors:
            all_valid = False
            for err in errors:
                console.print(f"  [red]✗[/red] AF3: {err}")
        else:
            console.print("  [green]✓[/green] AF3: all paths valid")

    if "boltz" in cfg.models_enabled:
        runner = BoltzRunner(cfg)
        errors = runner.validate_paths()
        if errors:
            all_valid = False
            for err in errors:
                console.print(f"  [red]✗[/red] Boltz: {err}")
        else:
            console.print("  [green]✓[/green] Boltz: all paths valid")

    if "rf3" in cfg.models_enabled:
        runner = RF3Runner(cfg)
        errors = runner.validate_paths()
        if errors:
            all_valid = False
            for err in errors:
                console.print(f"  [red]✗[/red] RF3: {err}")
        else:
            console.print("  [green]✓[/green] RF3: all paths valid")

    # Final status
    console.print()
    if all_valid:
        console.print("[green]All validations passed![/green]")
    else:
        console.print("[red]Some validations failed. Please fix the errors above.[/red]")
        raise typer.Exit(1)


@app.command()
def run(
    config: Annotated[
        Optional[Path],
        typer.Option("--config", "-c", help="Path to pipeline config YAML"),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", "-n", help="Show what would be submitted"),
    ] = False,
    resume: Annotated[
        bool,
        typer.Option("--resume", "-r", help="Skip completed jobs"),
    ] = True,
    model: Annotated[
        Optional[str],
        typer.Option("--model", "-m", help="Run only specific model (af3/boltz/rf3)"),
    ] = None,
    protein: Annotated[
        Optional[str],
        typer.Option("--protein", "-p", help="Run only specific protein ID"),
    ] = None,
    local: Annotated[
        bool,
        typer.Option("--local", help="Run locally instead of SLURM"),
    ] = False,
    # Model-specific overrides
    af3_seeds: Annotated[
        Optional[int],
        typer.Option("--af3-seeds", help="Override AF3 seed count"),
    ] = None,
    af3_diffusion_samples: Annotated[
        Optional[int],
        typer.Option("--af3-diffusion-samples", help="Override AF3 diffusion samples"),
    ] = None,
    boltz_diffusion_samples: Annotated[
        Optional[int],
        typer.Option("--boltz-diffusion-samples", help="Override Boltz diffusion samples"),
    ] = None,
    boltz_recycling_steps: Annotated[
        Optional[int],
        typer.Option("--boltz-recycling-steps", help="Override Boltz recycling steps"),
    ] = None,
    rf3_seed: Annotated[
        Optional[int],
        typer.Option("--rf3-seed", help="Override RF3 seed"),
    ] = None,
    rf3_diffusion_batch_size: Annotated[
        Optional[int],
        typer.Option("--rf3-diffusion-batch-size", help="Override RF3 diffusion batch size"),
    ] = None,
    rf3_n_recycles: Annotated[
        Optional[int],
        typer.Option("--rf3-n-recycles", help="Override RF3 recycle count"),
    ] = None,
    rf3_num_steps: Annotated[
        Optional[int],
        typer.Option("--rf3-num-steps", help="Override RF3 diffusion steps"),
    ] = None,
    ligand_ccd_list: Annotated[
        Optional[Path],
        typer.Option(
            "--ligand-ccd-list",
            help=(
                "Path to a .txt file with one CCD code per line. "
                "When provided, ligands are taken from this list and all "
                "models use standard CCD codes (no custom CIF paths). "
                "Must match the list used during 'manifest' generation."
            ),
        ),
    ] = None,
    oligo_definitions: Annotated[
        Optional[Path],
        typer.Option(
            "--oligo-definitions",
            help=(
                "Path to a YAML file defining oligosaccharide prefixes, "
                "monomer CCD codes, and bond atom pairs. "
                "When provided, AF3 builds multi-monomer ccdCodes and "
                "bondedAtomPairs for matching ligands instead of using "
                "custom CIF files. Overrides inputs.oligo_definitions "
                "in the config file."
            ),
        ),
    ] = None,
):
    """Run the structure prediction pipeline.

    Submits jobs grouped by (ligand, model).  For AF3 the two-stage
    fan-in is used: MSA per (protein, ligand) then inference per ligand
    with ``--dependency=afterok``.
    """
    cfg = get_config(config)

    # Apply CLI overrides
    if af3_seeds is not None:
        cfg.af3.seeds = af3_seeds
    if af3_diffusion_samples is not None:
        cfg.af3.num_diffusion_samples = af3_diffusion_samples
    if boltz_diffusion_samples is not None:
        cfg.boltz.diffusion_samples = boltz_diffusion_samples
    if boltz_recycling_steps is not None:
        cfg.boltz.recycling_steps = boltz_recycling_steps
    if rf3_seed is not None:
        cfg.rf3.seed = rf3_seed
    if rf3_diffusion_batch_size is not None:
        cfg.rf3.diffusion_batch_size = rf3_diffusion_batch_size
    if rf3_n_recycles is not None:
        cfg.rf3.n_recycles = rf3_n_recycles
    if rf3_num_steps is not None:
        cfg.rf3.num_steps = rf3_num_steps

    # ── Load manifests ─────────────────────────────────────────
    manifest_dir = cfg.outputs.manifest_dir
    cases_csv = manifest_dir / "cases.csv"
    proteins_csv = manifest_dir / "proteins.csv"
    msa_csv = manifest_dir / "msa.csv"
    ligands_csv = manifest_dir / "ligands.csv"

    if not cases_csv.exists():
        console.print("[red]Error:[/red] cases.csv not found. Run 'manifest' first.")
        raise typer.Exit(1)

    cases = load_cases_manifest(cases_csv)
    proteins = load_proteins_manifest(proteins_csv)
    protein_seqs = {p.protein_id: p.sequence for p in proteins}

    # MSA manifest
    msa_by_protein: dict[str, MSARecord] = {}
    if msa_csv.exists():
        msa_records = load_msa_manifest(msa_csv)
        msa_by_protein = {r.protein_id: r for r in msa_records}

    # Ligand CIF paths – only populated when CIF files exist (not CCD-list mode).
    # ligand_ccd_list flag controls whether we're in CCD-list mode.
    ligand_cif_paths: dict[str, Path] = {}
    if ligands_csv.exists() and ligand_ccd_list is None:
        ligands = load_ligands_manifest(ligands_csv)
        for lig in ligands:
            if lig.cif_path:  # non-empty means CIF file mode
                ligand_cif_paths[lig.ligand_id] = Path(lig.cif_path)

    # ── Filter ─────────────────────────────────────────────────
    if model:
        cases = [c for c in cases if c.model == model.lower()]
    if protein:
        cases = [c for c in cases if c.protein_id == protein]

    # ── Resume / pending ───────────────────────────────────────
    tracker = StatusTracker(cfg.outputs.work_dir)
    pending_cases = tracker.get_pending_cases(cases, resume=resume)

    if not pending_cases:
        console.print("[green]All cases are completed![/green]")
        return

    # ── Group by (ligand_ccd, model) ───────────────────────────
    grouped = group_cases_by_ligand_model(pending_cases)

    console.print(f"\n[bold]Jobs to submit:[/bold] {len(grouped)} groups")
    for (lig, mdl), group_cases in sorted(grouped.items()):
        pids = ", ".join(sorted(set(c.protein_id for c in group_cases)))
        console.print(f"  - {lig}/{mdl}: {len(group_cases)} proteins ({pids})")

    if dry_run:
        # Collect all paths that jobs would reference and show mount plan
        all_paths: list[str] = []
        all_paths.append(str(cfg.paths.af3_dir))
        all_paths.append(str(cfg.paths.af3_databases))
        all_paths.append(str(cfg.paths.boltz_image))
        all_paths.append(str(cfg.paths.boltz_weights))
        all_paths.append(str(cfg.paths.rf3_foundry_root))
        all_paths.append(str(cfg.paths.rf3_image))
        all_paths.append(str(cfg.paths.rf3_checkpoint))
        if cfg.inputs.msa_dir:
            all_paths.append(str(cfg.inputs.msa_dir))
        for sqsh in cfg.inputs.msa_sqsh_files:
            all_paths.append(str(sqsh))
        for cif_path in ligand_cif_paths.values():
            all_paths.append(str(cif_path))
        all_paths.append(str(cfg.outputs.work_dir))

        mounts = compute_bind_mounts(all_paths)
        console.print(f"\n{log_mount_plan(mounts, all_paths)}")
        console.print("\n[yellow]Dry run mode — no jobs submitted[/yellow]")
        return

    # ── Load oligo registry ────────────────────────────────
    oligo_registry = OligoRegistry.empty()
    oligo_path = oligo_definitions or cfg.inputs.oligo_definitions
    if oligo_path is not None:
        oligo_path = Path(oligo_path)
        if not oligo_path.is_absolute():
            oligo_path = oligo_path.resolve()
        if not oligo_path.exists():
            console.print(
                f"[red]Error:[/red] oligo definitions file not found: {oligo_path}"
            )
            raise typer.Exit(1)
        try:
            oligo_registry = OligoRegistry.from_yaml(oligo_path)
            console.print(
                f"  Loaded oligo definitions from {oligo_path} "
                f"({len(oligo_registry.specs)} prefixes)"
            )
        except Exception as e:
            console.print(f"[red]Error loading oligo definitions:[/red] {e}")
            raise typer.Exit(1)

    # ── Setup ──────────────────────────────────────────────────
    executor = LocalExecutor(dry_run=dry_run) if local else SlurmExecutor(dry_run=dry_run)
    runners = {
        "af3": AF3Runner(cfg, oligo_registry=oligo_registry),
        "boltz": BoltzRunner(cfg),
        "rf3": RF3Runner(cfg),
    }

    console.print("\n[bold]Submitting jobs...[/bold]")
    submitted_jobs = []

    # ── Collect AF3 MSA jobs first (needed for dependency) ─────
    # One MSA job per protein; reuse across ligands.
    af3_msa_jobs: dict[str, dict] = {}  # key = protein_id
    af3_msa_ligand: dict[str, str] = {}

    for (ligand_ccd, model_name), group_cases in sorted(grouped.items()):
        if model_name != "af3":
            continue
        for case in group_cases:
            af3_msa_ligand.setdefault(case.protein_id, ligand_ccd)

    for (ligand_ccd, model_name), group_cases in sorted(grouped.items()):
        if model_name != "af3":
            continue

        af3_runner: AF3Runner = runners["af3"]  # type: ignore[assignment]

        for case in group_cases:
            msa_key = case.protein_id
            if msa_key in af3_msa_jobs:
                continue  # already submitted

            msa_work_dir = cfg.outputs.work_dir / "af3_msa" / msa_key
            msa_job_name = f"af3_msa_{msa_key}"
            msa_ligand = af3_msa_ligand.get(case.protein_id, ligand_ccd)

            # Get CIF path for this ligand
            ligand_cif_path = ligand_cif_paths.get(case.ligand_id)

            msa_script = af3_runner.build_msa_slurm_script(
                protein_id=case.protein_id,
                protein_sequence=protein_seqs[case.protein_id],
                ligand_ccd_code=msa_ligand,
                work_dir=msa_work_dir,
                job_name=msa_job_name,
                ligand_cif_path=ligand_cif_path,
            )

            msa_job = executor.submit(msa_script, msa_work_dir, msa_job_name)
            submitted_jobs.append(msa_job)
            af3_msa_jobs[msa_key] = {
                "job": msa_job,
                "work_dir": msa_work_dir,
            }
            console.print(f"  [green]✓[/green] {msa_job.job_name} → {msa_job.job_id}")

    # ── Submit (ligand, model) jobs ────────────────────────────
    for (ligand_ccd, model_name), group_cases in sorted(grouped.items()):
        runner = runners[model_name]
        work_dir = cfg.outputs.work_dir / ligand_ccd / model_name
        job_name = f"{model_name}_{ligand_ccd}"

        # Build protein_sequences dict for this group
        group_protein_seqs = {
            c.protein_id: protein_seqs[c.protein_id] for c in group_cases
        }

        if model_name == "af3":
            # ── AF3 inference (depends on MSA jobs) ────────────
            af3_runner = runners["af3"]  # type: ignore

            # Collect MSA work dirs and dependency job IDs
            msa_work_dirs: dict[str, Path] = {}
            dep_job_ids: list[str] = []

            for case in group_cases:
                msa_key = case.protein_id
                msa_info = af3_msa_jobs.get(msa_key)
                if msa_info:
                    msa_work_dirs[case.protein_id] = msa_info["work_dir"]
                    dep_job_ids.append(msa_info["job"].job_id)

            # Build inference script
            inf_script = af3_runner.build_slurm_script(
                cases=group_cases,
                protein_sequences=group_protein_seqs,
                work_dir=work_dir,
                job_name=f"{job_name}_inf",
                msa_work_dirs=msa_work_dirs,
                ligand_cif_paths=ligand_cif_paths,
            )

            # Submit with dependency on all MSA jobs for this ligand
            dep_str = None
            if dep_job_ids and not local:
                dep_str = "afterok:" + ":".join(dep_job_ids)

            inf_job = executor.submit(inf_script, work_dir, f"{job_name}_inf", dep_str)
            submitted_jobs.append(inf_job)
            console.print(f"  [green]✓[/green] {inf_job.job_name} → {inf_job.job_id}")

        elif model_name == "rf3":
            # ── RF3 (single job per ligand) ────────────────────
            # Collect MSA paths for all proteins
            msa_paths: dict[str, str] = {}
            msa_sqsh_file = None
            for case in group_cases:
                rec = msa_by_protein.get(case.protein_id)
                if rec:
                    if rec.sqsh_file:
                        msa_sqsh_file = Path(rec.sqsh_file)
                    elif rec.msa_path:
                        msa_paths[case.protein_id] = rec.msa_path

            script = runner.build_slurm_script(
                cases=group_cases,
                protein_sequences=group_protein_seqs,
                work_dir=work_dir,
                job_name=job_name,
                msa_paths=msa_paths,
                msa_sqsh_file=msa_sqsh_file,
                ligand_cif_paths=ligand_cif_paths,
            )
            job = executor.submit(script, work_dir, job_name)
            submitted_jobs.append(job)
            console.print(f"  [green]✓[/green] {job.job_name} → {job.job_id}")

        else:
            # ── Boltz (single job per ligand) ──────────────────
            msa_paths_boltz: dict[str, str] = {}
            msa_sqsh_file_boltz = None
            for case in group_cases:
                rec = msa_by_protein.get(case.protein_id)
                if rec:
                    if rec.sqsh_file:
                        msa_sqsh_file_boltz = Path(rec.sqsh_file)
                    elif rec.msa_path:
                        msa_paths_boltz[case.protein_id] = rec.msa_path

            script = runner.build_slurm_script(
                cases=group_cases,
                protein_sequences=group_protein_seqs,
                work_dir=work_dir,
                job_name=job_name,
                msa_paths=msa_paths_boltz,
                msa_sqsh_file=msa_sqsh_file_boltz,
            )
            job = executor.submit(script, work_dir, job_name)
            submitted_jobs.append(job)
            console.print(f"  [green]✓[/green] {job.job_name} → {job.job_id}")

    console.print(f"\n[bold]Submitted {len(submitted_jobs)} jobs[/bold]")


@app.command()
def status(
    config: Annotated[
        Optional[Path],
        typer.Option("--config", "-c", help="Path to pipeline config YAML"),
    ] = None,
):
    """Show status of all jobs."""
    cfg = get_config(config)
    tracker = StatusTracker(cfg.outputs.work_dir)

    summary = tracker.generate_summary()

    console.print(f"\n[bold]Pipeline Status[/bold]")
    console.print(f"Total ligands: {summary['total_ligands']}")

    if summary["by_model"]:
        table = Table()
        table.add_column("Model")
        table.add_column("Completed", justify="right")
        table.add_column("Pending", justify="right")
        table.add_column("Failed", justify="right")

        for model, counts in sorted(summary["by_model"].items()):
            table.add_row(
                model,
                str(counts.get("completed", 0)),
                str(counts.get("pending", 0)),
                str(counts.get("failed", 0)),
            )

        console.print(table)
    else:
        console.print("[yellow]No jobs have been run yet[/yellow]")


@app.command()
def init_config(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output path for config file"),
    ] = Path("config/pipeline.yaml"),
):
    """Generate an example configuration file."""
    example_config = '''# Structure Pipeline Configuration
# See README.md for full documentation

paths:
  # AF3 paths
  af3_dir: /cluster/projects/nn1003k/prog/af3
  af3_cpu_image: /cluster/projects/nn1003k/prog/af3/af3_cpu_amd64.sif
  af3_gpu_image: /cluster/projects/nn1003k/prog/af3/af3_gpu_arm64.sif
  af3_weights: /cluster/projects/nn1003k/prog/af3/weights
  af3_databases: /cluster/work/shared/alphafold_uncompressed.squashfs

  # Boltz paths
  boltz_image: /cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif
  boltz_weights: /cluster/projects/nn1003k/prog/boltz/weights

  # RF3 paths
  rf3_foundry_root: /cluster/projects/nn1003k/prog/foundry
  rf3_image: /cluster/projects/nn1003k/prog/foundry/foundry.sif
  rf3_checkpoint: /cluster/projects/nn1003k/prog/foundry/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt

slurm:
  account: nn1003k
  partition_cpu: normal
  partition_gpu: accel
    mem_per_gpu: 80G
    # Per-model memory overrides (uncomment to customise):
    # mem_per_gpu_af3: 80G
    # mem_per_gpu_boltz: 80G
    # mem_per_gpu_rf3: 80G
  time_msa: "02:00:00"
  time_inference: "01:00:00"
  gpus: 1
  cpus_msa: 8
  mem_per_cpu_msa: 10G
  tmp_base: /cluster/work/projects/nn1003k/eirik/tmp

af3:
  seeds: 10
  num_diffusion_samples: 5

boltz:
  seeds: 5
  diffusion_samples: 5
  recycling_steps: 10
  use_affinity: true
  use_msa_server: false

rf3:
  seed: 42
  diffusion_batch_size: 5
  n_recycles: 10
  num_steps: 50
  early_stopping_plddt_threshold: 0.5
  skip_existing: false

inputs:
  proteins_fasta: input/fasta/proteins.fasta
  ligands_dir: input/ligand
  msa_dir: input/msa
  msa_sqsh_files: []
  msa_sqsh_mount_point: /msa

outputs:
  work_dir: work
  manifest_dir: manifest
  results_dir: results

models_enabled:
  - af3
  - boltz
  - rf3
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as f:
        f.write(example_config)

    console.print(f"[green]✓[/green] Example config written to {output}")
    console.print("\nEdit the paths section to match your environment.")


def main():
    """Entry point."""
    app()


if __name__ == "__main__":
    main()
