from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lpmo_pipeline.analysis.convergence_metrics import (
    ConditionConvergenceSummary,
    ConvergenceConfig,
    ConvergencePoseInput,
    PoseConvergenceMetrics,
    build_condition_convergence_summary,
    compute_condition_convergence,
    select_reference_pose,
    write_condition_convergence_summary_tsv,
    write_pose_convergence_tsv,
)


def _pdb_atom_line(
    record: str,
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.00:>6.2f}{20.00:>6.2f}          {element:>2}"
    )


def _write_test_complex(
    output_path: Path,
    *,
    protein_shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ligand_shift: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Path:
    protein_coords = [
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (0.0, 1.5, 0.0),
    ]
    ligand_coords = [
        (2.0, 2.0, 0.0),
        (3.0, 2.0, 0.0),
    ]
    px, py, pz = protein_shift
    lx, ly, lz = ligand_shift

    lines = [
        _pdb_atom_line("ATOM", 1, "CA", "ALA", "A", 1, protein_coords[0][0] + px, protein_coords[0][1] + py, protein_coords[0][2] + pz, "C"),
        _pdb_atom_line("ATOM", 2, "CA", "GLY", "A", 2, protein_coords[1][0] + px, protein_coords[1][1] + py, protein_coords[1][2] + pz, "C"),
        _pdb_atom_line("ATOM", 3, "CA", "SER", "A", 3, protein_coords[2][0] + px, protein_coords[2][1] + py, protein_coords[2][2] + pz, "C"),
        _pdb_atom_line("HETATM", 4, "C1", "NAG", "B", 1, ligand_coords[0][0] + px + lx, ligand_coords[0][1] + py + ly, ligand_coords[0][2] + pz + lz, "C"),
        _pdb_atom_line("HETATM", 5, "O4", "NAG", "B", 1, ligand_coords[1][0] + px + lx, ligand_coords[1][1] + py + ly, ligand_coords[1][2] + pz + lz, "O"),
        "END",
    ]
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def test_select_reference_pose_prefers_seed_one_candidates() -> None:
    pose_inputs = [
        ConvergencePoseInput(
            pose_id="pose-seed2",
            condition_id="cond",
            complex_pdb=Path("/tmp/pose-seed2.pdb"),
            seed=2,
            sample=0,
        ),
        ConvergencePoseInput(
            pose_id="pose-seed1",
            condition_id="cond",
            complex_pdb=Path("/tmp/pose-seed1.pdb"),
            seed=1,
            sample=3,
        ),
    ]

    reference = select_reference_pose(pose_inputs)

    assert reference.pose_id == "pose-seed1"


def test_select_reference_pose_prefers_higher_ranking_within_seed_one() -> None:
    pose_inputs = [
        ConvergencePoseInput(
            pose_id="pose-low-rank",
            condition_id="cond",
            complex_pdb=Path("/tmp/pose-low-rank.pdb"),
            seed=1,
            sample=0,
            ranking_score=0.70,
        ),
        ConvergencePoseInput(
            pose_id="pose-high-rank",
            condition_id="cond",
            complex_pdb=Path("/tmp/pose-high-rank.pdb"),
            seed=1,
            sample=5,
            ranking_score=0.95,
        ),
    ]

    reference = select_reference_pose(pose_inputs)

    assert reference.pose_id == "pose-high-rank"


def test_compute_condition_convergence_aligns_protein_before_ligand_rmsd(tmp_path: Path) -> None:
    ref_path = _write_test_complex(tmp_path / "ref.pdb")
    shifted_path = _write_test_complex(
        tmp_path / "shifted.pdb",
        protein_shift=(10.0, 4.0, 1.0),
        ligand_shift=(0.0, 0.0, 0.0),
    )

    pose_metrics, summary = compute_condition_convergence(
        [
            ConvergencePoseInput(
                pose_id="pose-ref",
                condition_id="cond",
                complex_pdb=ref_path,
                seed=1,
                sample=0,
            ),
            ConvergencePoseInput(
                pose_id="pose-shifted",
                condition_id="cond",
                complex_pdb=shifted_path,
                seed=2,
                sample=0,
            ),
        ],
        config=ConvergenceConfig(convergent_rmsd_max_a=0.5),
    )

    assert [metric.reference_pose_id for metric in pose_metrics] == ["pose-ref", "pose-ref"]
    assert pose_metrics[0].ligand_rmsd_to_reference == pytest.approx(0.0, abs=1e-6)
    assert pose_metrics[1].ligand_rmsd_to_reference == pytest.approx(0.0, abs=1e-4)
    assert pose_metrics[1].convergent_flag is True
    assert summary.convergence_fraction == pytest.approx(1.0)


def test_compute_condition_convergence_detects_shifted_ligand(tmp_path: Path) -> None:
    ref_path = _write_test_complex(tmp_path / "ref.pdb")
    shifted_ligand_path = _write_test_complex(
        tmp_path / "shifted_ligand.pdb",
        protein_shift=(5.0, 0.0, 0.0),
        ligand_shift=(3.0, 0.0, 0.0),
    )

    pose_metrics, summary = compute_condition_convergence(
        [
            ConvergencePoseInput(
                pose_id="pose-ref",
                condition_id="cond",
                complex_pdb=ref_path,
                seed=1,
                sample=0,
            ),
            ConvergencePoseInput(
                pose_id="pose-shifted",
                condition_id="cond",
                complex_pdb=shifted_ligand_path,
                seed=2,
                sample=0,
            ),
        ],
        config=ConvergenceConfig(convergent_rmsd_max_a=2.0),
    )

    assert pose_metrics[1].ligand_rmsd_to_reference == pytest.approx(3.0, abs=1e-4)
    assert pose_metrics[1].convergent_flag is False
    assert summary.convergence_fraction == pytest.approx(0.5)
    assert summary.low_convergence_flag is False


def test_build_condition_convergence_summary_aggregates_iqr() -> None:
    summary = build_condition_convergence_summary(
        [
            PoseConvergenceMetrics("pose-1", "cond", "pose-1", 0.0, True),
            PoseConvergenceMetrics("pose-2", "cond", "pose-1", 1.0, True),
            PoseConvergenceMetrics("pose-3", "cond", "pose-1", 3.0, False),
        ],
        config=ConvergenceConfig(low_convergence_flag_threshold=0.8),
    )

    assert summary == ConditionConvergenceSummary(
        condition_id="cond",
        reference_pose_id="pose-1",
        n_qc_pass_poses=3,
        convergence_fraction=pytest.approx(2 / 3),
        median_ligand_rmsd=pytest.approx(1.0),
        iqr_ligand_rmsd=pytest.approx(1.5),
        low_convergence_flag=True,
    )


def test_convergence_tsv_writers_emit_expected_columns(tmp_path: Path) -> None:
    pose_path = tmp_path / "pose_convergence.tsv"
    condition_path = tmp_path / "condition_convergence.tsv"
    write_pose_convergence_tsv(
        [PoseConvergenceMetrics("pose-1", "cond", "pose-1", 0.25, True)],
        pose_path,
    )
    write_condition_convergence_summary_tsv(
        [ConditionConvergenceSummary("cond", "pose-1", 1, 1.0, 0.25, 0.0, False)],
        condition_path,
    )

    with pose_path.open(newline="") as handle:
        pose_rows = list(csv.DictReader(handle, delimiter="\t"))
    with condition_path.open(newline="") as handle:
        condition_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert pose_rows == [
        {
            "pose_id": "pose-1",
            "condition_id": "cond",
            "reference_pose_id": "pose-1",
            "ligand_rmsd_to_reference": "0.25",
            "convergent_flag": "True",
        }
    ]
    assert condition_rows == [
        {
            "condition_id": "cond",
            "reference_pose_id": "pose-1",
            "n_qc_pass_poses": "1",
            "convergence_fraction": "1.0",
            "median_ligand_rmsd": "0.25",
            "iqr_ligand_rmsd": "0.0",
            "low_convergence_flag": "False",
        }
    ]