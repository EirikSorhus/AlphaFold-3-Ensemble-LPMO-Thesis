# tests/test_discovery.py
"""
Tests for lpmo_pipeline.io.discovery – work-root artifact discovery.

Covers:
  - Target discovery by naming format (CEL6, STA8, NAG4, …)
  - AF3/RF3 inclusion, Boltz exclusion
  - Seed/sample parsing
  - Empty directories
  - Missing CIF, missing JSON
  - Discovery of all seeds/samples
  - Complete vs incomplete entries
  - Multiple runs per model
  - Uniprot-target directory parsing
  - WorkRootManifest iteration helpers
"""

from __future__ import annotations

from pathlib import Path

import pytest

from lpmo_pipeline.io.discovery import (
    TARGET_RE,
    SEED_SAMPLE_RE,
    UNIPROT_TARGET_RE,
    RUN_ID_RE,
    discover_work_root,
    parse_seed_sample,
)
from lpmo_pipeline.utils.data_models import (
    ALLOWED_MODELS,
    SampleStatus,
    WorkRootManifest,
)


# ---------------------------------------------------------------------------
# Fixtures: build minimal work-root trees in tmp_path
# ---------------------------------------------------------------------------

def _touch(path: Path) -> Path:
    """Create a file and its parents."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


def _make_sample_dir(parent: Path, seed: int, sample: int) -> Path:
    """Create empty seed-sample directory."""
    d = parent / f"seed-{seed}_sample-{sample}"
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture
def minimal_work(tmp_path: Path) -> Path:
    """Minimal valid work root: one target, one model, one run, one uniprot."""
    root = tmp_path / "work"
    ut = root / "CEL6" / "af3" / "runs" / "100001" / "B6EQJ6_CEL6"
    _touch(ut / "B6EQJ6_CEL6_model.cif")
    _make_sample_dir(ut, seed=1, sample=0)
    _make_sample_dir(ut, seed=1, sample=1)
    return root


@pytest.fixture
def full_work(tmp_path: Path) -> Path:
    """Full work root with multiple targets, models, runs, seed/samples."""
    root = tmp_path / "work"

    # CEL6 — af3 + rf3
    for model, run_id in [("af3", "187030"), ("rf3", "182585")]:
        for uid in ["B6EQJ6", "D0EW65"]:
            ut = root / "CEL6" / model / "runs" / run_id / f"{uid}_CEL6"
            _touch(ut / f"{uid}_CEL6_model.cif")
            for s in range(1, 3):
                for m in range(2):
                    _make_sample_dir(ut, seed=s, sample=m)

    # STA6 — af3 only (no CIF yet, just seed/sample dirs)
    ut = root / "STA6" / "af3" / "runs" / "187032" / "D5UGB1_STA6"
    for s in range(1, 4):
        for m in range(3):
            _make_sample_dir(ut, seed=s, sample=m)

    # NAG6 — rf3 with a CIF + confidence JSON in seed/sample dir
    ut = root / "NAG6" / "rf3" / "runs" / "182587" / "Q9RFX5_NAG6"
    _touch(ut / "Q9RFX5_NAG6_model.cif")
    sd = _make_sample_dir(ut, seed=42, sample=0)
    _touch(sd / "Q9RFX5_NAG6_seed42_sample0.cif")
    _touch(sd / "Q9RFX5_NAG6_seed42_sample0_confidence.json")

    # Boltz should be ignored
    bt = root / "CEL6" / "boltz" / "runs" / "999" / "B6EQJ6_CEL6"
    _touch(bt / "B6EQJ6_CEL6_model.cif")

    # Non-target directory should be ignored
    (root / "af3_msa").mkdir(parents=True, exist_ok=True)

    return root


# ---------------------------------------------------------------------------
# Regex unit tests
# ---------------------------------------------------------------------------

class TestRegexPatterns:
    """Validate the regex patterns used in discovery."""

    @pytest.mark.parametrize("name", ["CEL6", "STA8", "NAG4", "AMY12", "XY1234"])
    def test_target_re_valid(self, name: str) -> None:
        assert TARGET_RE.match(name) is not None

    @pytest.mark.parametrize("name", ["af3_msa", "cache", "boltz", "cel6", "123", "", "A", "AB"])
    def test_target_re_invalid(self, name: str) -> None:
        assert TARGET_RE.match(name) is None

    @pytest.mark.parametrize("name,seed,sample", [
        ("seed-1_sample-0", 1, 0),
        ("seed-10_sample-4", 10, 4),
        ("seed-42_sample-0", 42, 0),
        ("seed-100_sample-99", 100, 99),
    ])
    def test_seed_sample_re_valid(self, name: str, seed: int, sample: int) -> None:
        m = SEED_SAMPLE_RE.match(name)
        assert m is not None
        assert int(m.group(1)) == seed
        assert int(m.group(2)) == sample

    @pytest.mark.parametrize("name", ["seed-_sample-0", "sample-0_seed-1", "cache", ""])
    def test_seed_sample_re_invalid(self, name: str) -> None:
        assert SEED_SAMPLE_RE.match(name) is None

    @pytest.mark.parametrize("name,uid,target", [
        ("B6EQJ6_CEL6", "B6EQJ6", "CEL6"),
        ("D0EW65_STA8", "D0EW65", "STA8"),
        ("A0A1B2C3D4_NAG12", "A0A1B2C3D4", "NAG12"),
    ])
    def test_uniprot_target_re_valid(self, name: str, uid: str, target: str) -> None:
        m = UNIPROT_TARGET_RE.match(name)
        assert m is not None
        assert m.group(1) == uid
        assert m.group(2) == target

    @pytest.mark.parametrize("name", ["cache", "short_C6", "", "B6_CEL"])
    def test_uniprot_target_re_invalid(self, name: str) -> None:
        assert UNIPROT_TARGET_RE.match(name) is None

    @pytest.mark.parametrize("name", ["187030", "1", "9999999"])
    def test_run_id_re_valid(self, name: str) -> None:
        assert RUN_ID_RE.match(name) is not None

    @pytest.mark.parametrize("name", ["abc", "run-1", ""])
    def test_run_id_re_invalid(self, name: str) -> None:
        assert RUN_ID_RE.match(name) is None


class TestParseSeedSample:
    def test_valid(self) -> None:
        assert parse_seed_sample("seed-3_sample-2") == (3, 2)

    def test_invalid(self) -> None:
        assert parse_seed_sample("not-a-seed") is None


# ---------------------------------------------------------------------------
# Allowed models constant
# ---------------------------------------------------------------------------

class TestAllowedModels:
    def test_af3_included(self) -> None:
        assert "af3" in ALLOWED_MODELS

    def test_rf3_included(self) -> None:
        assert "rf3" in ALLOWED_MODELS

    def test_boltz_excluded(self) -> None:
        assert "boltz" not in ALLOWED_MODELS
        assert "boltz2" not in ALLOWED_MODELS


# ---------------------------------------------------------------------------
# Discovery: minimal work root
# ---------------------------------------------------------------------------

class TestDiscoverMinimal:
    def test_discovers_one_target(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        assert len(m.targets) == 1
        assert m.targets[0].target == "CEL6"

    def test_discovers_af3_model(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        models = [mdl.model for mdl in m.targets[0].models]
        assert models == ["af3"]

    def test_discovers_run_id(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        runs = m.targets[0].models[0].runs
        assert len(runs) == 1
        assert runs[0].run_id == "100001"

    def test_discovers_uniprot_target(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        ut = m.targets[0].models[0].runs[0].uniprot_targets[0]
        assert ut.uniprot_id == "B6EQJ6"
        assert ut.target == "CEL6"

    def test_model_cif_found(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        ut = m.targets[0].models[0].runs[0].uniprot_targets[0]
        assert ut.has_model_cif
        assert ut.model_cif_path.name == "B6EQJ6_CEL6_model.cif"

    def test_seed_sample_dirs_discovered(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        ut = m.targets[0].models[0].runs[0].uniprot_targets[0]
        assert ut.total_samples == 2
        assert ut.seeds == [1]
        assert ut.sample_indices == [0, 1]

    def test_empty_seed_sample_status(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        ut = m.targets[0].models[0].runs[0].uniprot_targets[0]
        for s in ut.samples:
            assert s.status == SampleStatus.EMPTY_DIRECTORY

    def test_no_errors(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        assert m.errors == []


# ---------------------------------------------------------------------------
# Discovery: full work root
# ---------------------------------------------------------------------------

class TestDiscoverFull:
    def test_three_targets_discovered(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        target_names = sorted(t.target for t in m.targets)
        assert target_names == ["CEL6", "NAG6", "STA6"]

    def test_boltz_excluded(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        cel6 = [t for t in m.targets if t.target == "CEL6"][0]
        model_names = [mdl.model for mdl in cel6.models]
        assert "boltz" not in model_names

    def test_af3_and_rf3_for_cel6(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        cel6 = [t for t in m.targets if t.target == "CEL6"][0]
        model_names = sorted(mdl.model for mdl in cel6.models)
        assert model_names == ["af3", "rf3"]

    def test_af3_msa_dir_ignored(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        target_names = [t.target for t in m.targets]
        assert "af3_msa" not in target_names

    def test_cel6_af3_uniprots(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        cel6 = [t for t in m.targets if t.target == "CEL6"][0]
        af3 = [mdl for mdl in cel6.models if mdl.model == "af3"][0]
        ut_ids = sorted(ut.uniprot_id for r in af3.runs for ut in r.uniprot_targets)
        assert ut_ids == ["B6EQJ6", "D0EW65"]

    def test_sta6_no_model_cif(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        sta6 = [t for t in m.targets if t.target == "STA6"][0]
        af3 = sta6.models[0]
        ut = af3.runs[0].uniprot_targets[0]
        assert not ut.has_model_cif
        assert ut.model_cif_path is None

    def test_sta6_seeds_and_samples(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        sta6 = [t for t in m.targets if t.target == "STA6"][0]
        ut = sta6.models[0].runs[0].uniprot_targets[0]
        assert ut.seeds == [1, 2, 3]
        assert ut.sample_indices == [0, 1, 2]
        assert ut.total_samples == 9  # 3 seeds × 3 samples

    def test_nag6_complete_sample(self, full_work: Path) -> None:
        """Sample with both CIF and JSON → COMPLETE status."""
        m = discover_work_root(full_work)
        nag6 = [t for t in m.targets if t.target == "NAG6"][0]
        ut = nag6.models[0].runs[0].uniprot_targets[0]
        complete = [s for s in ut.samples if s.status == SampleStatus.COMPLETE]
        assert len(complete) == 1
        assert complete[0].seed == 42
        assert complete[0].sample == 0
        assert len(complete[0].cif_paths) == 1
        assert len(complete[0].confidence_json_paths) == 1

    def test_iter_all_model_cifs(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        cifs = list(m.iter_all_model_cifs())
        # CEL6 af3: 2 uniprots, CEL6 rf3: 2, NAG6 rf3: 1 top-level + 1 sample
        assert len(cifs) >= 5

    def test_summary_counts(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        s = m.summary
        assert s["targets"] == 3
        assert s["total_cif_files"] >= 5
        assert s["errors"] == 0


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_nonexistent_work_root(self, tmp_path: Path) -> None:
        m = discover_work_root(tmp_path / "does_not_exist")
        assert len(m.targets) == 0
        assert len(m.errors) == 1
        assert "not a directory" in m.errors[0]

    def test_empty_work_root(self, tmp_path: Path) -> None:
        root = tmp_path / "work"
        root.mkdir()
        m = discover_work_root(root)
        assert len(m.targets) == 0
        assert m.errors == []

    def test_target_with_no_models(self, tmp_path: Path) -> None:
        """Target directory exists but has no af3/rf3 subdirectories."""
        root = tmp_path / "work"
        (root / "CEL6").mkdir(parents=True)
        m = discover_work_root(root)
        assert len(m.targets) == 0  # excluded because no models

    def test_model_with_no_runs_dir(self, tmp_path: Path) -> None:
        """af3 directory exists but has no runs/ subdirectory."""
        root = tmp_path / "work"
        (root / "CEL6" / "af3").mkdir(parents=True)
        m = discover_work_root(root)
        cel6 = m.targets[0]
        assert len(cel6.models) == 1
        assert len(cel6.models[0].runs) == 0

    def test_run_with_no_uniprot_dirs(self, tmp_path: Path) -> None:
        """Run directory exists but is empty."""
        root = tmp_path / "work"
        (root / "CEL6" / "af3" / "runs" / "999").mkdir(parents=True)
        m = discover_work_root(root)
        run = m.targets[0].models[0].runs[0]
        assert len(run.uniprot_targets) == 0

    def test_cache_dir_ignored_in_run(self, tmp_path: Path) -> None:
        """RF3 'cache' directory inside a run should be ignored."""
        root = tmp_path / "work"
        run_dir = root / "CEL6" / "rf3" / "runs" / "100"
        (run_dir / "cache").mkdir(parents=True)
        ut = run_dir / "B6EQJ6_CEL6"
        _touch(ut / "model.cif")
        m = discover_work_root(root)
        run = m.targets[0].models[0].runs[0]
        assert len(run.uniprot_targets) == 1
        assert run.uniprot_targets[0].uniprot_id == "B6EQJ6"

    def test_missing_cif_only_json_in_sample(self, tmp_path: Path) -> None:
        """Sample directory with JSON but no CIF → MISSING_CIF."""
        root = tmp_path / "work"
        ut = root / "CEL6" / "af3" / "runs" / "100" / "B6EQJ6_CEL6"
        sd = _make_sample_dir(ut, seed=1, sample=0)
        _touch(sd / "confidence.json")
        m = discover_work_root(root)
        sample = m.targets[0].models[0].runs[0].uniprot_targets[0].samples[0]
        assert sample.status == SampleStatus.MISSING_CIF

    def test_cif_without_json_in_sample(self, tmp_path: Path) -> None:
        """Sample directory with CIF but no JSON → MISSING_CONFIDENCE_JSON."""
        root = tmp_path / "work"
        ut = root / "CEL6" / "af3" / "runs" / "100" / "B6EQJ6_CEL6"
        sd = _make_sample_dir(ut, seed=1, sample=0)
        _touch(sd / "output.cif")
        m = discover_work_root(root)
        sample = m.targets[0].models[0].runs[0].uniprot_targets[0].samples[0]
        assert sample.status == SampleStatus.MISSING_CONFIDENCE_JSON

    def test_multiple_runs_same_model(self, tmp_path: Path) -> None:
        """Multiple numeric run dirs under one model."""
        root = tmp_path / "work"
        for rid in ["100", "200", "300"]:
            ut = root / "CEL6" / "af3" / "runs" / rid / "B6EQJ6_CEL6"
            _touch(ut / "model.cif")
        m = discover_work_root(root)
        runs = m.targets[0].models[0].runs
        assert len(runs) == 3
        run_ids = [r.run_id for r in runs]
        assert sorted(run_ids) == ["100", "200", "300"]


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

class TestSerialization:
    def test_to_json_roundtrip(self, minimal_work: Path) -> None:
        """Manifest serializes to valid JSON."""
        import json
        m = discover_work_root(minimal_work)
        text = m.to_json()
        parsed = json.loads(text)
        assert parsed["summary"]["targets"] == 1

    def test_to_dict_has_all_keys(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        d = m.to_dict()
        assert "work_root" in d
        assert "summary" in d
        assert "errors" in d
        assert "targets" in d

    def test_sample_entry_to_dict(self, minimal_work: Path) -> None:
        m = discover_work_root(minimal_work)
        ut = m.targets[0].models[0].runs[0].uniprot_targets[0]
        sd = ut.samples[0].to_dict()
        assert "seed" in sd
        assert "sample" in sd
        assert "status" in sd
        assert sd["status"] == "empty_directory"


# ---------------------------------------------------------------------------
# Filtering: af3_only and latest_only
# ---------------------------------------------------------------------------

class TestAf3OnlyFilter:
    """Tests for the af3_only parameter."""

    def test_af3_only_skips_rf3(self, full_work: Path) -> None:
        m = discover_work_root(full_work, af3_only=True)
        all_models = [
            mdl.model for t in m.targets for mdl in t.models
        ]
        assert "rf3" not in all_models
        assert "af3" in all_models

    def test_af3_only_excludes_rf3_only_target(self, full_work: Path) -> None:
        """NAG6 only has rf3 in full_work → should be excluded entirely."""
        m = discover_work_root(full_work, af3_only=True)
        target_names = [t.target for t in m.targets]
        assert "NAG6" not in target_names

    def test_af3_only_keeps_af3_targets(self, full_work: Path) -> None:
        m = discover_work_root(full_work, af3_only=True)
        target_names = sorted(t.target for t in m.targets)
        assert "CEL6" in target_names
        assert "STA6" in target_names

    def test_default_includes_rf3(self, full_work: Path) -> None:
        m = discover_work_root(full_work)
        all_models = [mdl.model for t in m.targets for mdl in t.models]
        assert "rf3" in all_models


class TestLatestOnlyFilter:
    """Tests for the latest_only parameter."""

    @pytest.fixture
    def work_with_latest(self, tmp_path: Path) -> Path:
        """Work root with multiple runs and a 'latest' symlink."""
        root = tmp_path / "work"
        af3_dir = root / "CEL6" / "af3"
        runs_dir = af3_dir / "runs"

        # Two runs: old and new
        for rid, uid in [("100", "B6EQJ6"), ("200", "B6EQJ6")]:
            ut = runs_dir / rid / f"{uid}_CEL6"
            _touch(ut / f"{uid}_CEL6_model.cif")
            sd = _make_sample_dir(ut, seed=1, sample=0)
            _touch(sd / "pose.cif")
            _touch(sd / "confidence.json")

        # latest → 200
        latest = af3_dir / "latest"
        latest.symlink_to(runs_dir / "200")

        return root

    def test_latest_only_single_run(self, work_with_latest: Path) -> None:
        m = discover_work_root(work_with_latest, latest_only=True)
        runs = m.targets[0].models[0].runs
        assert len(runs) == 1
        assert runs[0].run_id == "200"

    def test_without_latest_flag_all_runs(self, work_with_latest: Path) -> None:
        m = discover_work_root(work_with_latest, latest_only=False)
        runs = m.targets[0].models[0].runs
        assert len(runs) == 2

    def test_latest_only_no_symlink_falls_back(self, tmp_path: Path) -> None:
        """Without a 'latest' symlink, latest_only scans all runs."""
        root = tmp_path / "work"
        for rid in ["100", "200"]:
            ut = root / "CEL6" / "af3" / "runs" / rid / "B6EQJ6_CEL6"
            _touch(ut / "model.cif")
        m = discover_work_root(root, latest_only=True)
        runs = m.targets[0].models[0].runs
        assert len(runs) == 2

    def test_combined_af3_only_latest_only(self, work_with_latest: Path) -> None:
        """Both flags together: only AF3, only latest run."""
        # Add an rf3 model to verify it's excluded
        rf3_ut = work_with_latest / "CEL6" / "rf3" / "runs" / "300" / "B6EQJ6_CEL6"
        _touch(rf3_ut / "model.cif")

        m = discover_work_root(
            work_with_latest, af3_only=True, latest_only=True,
        )
        assert len(m.targets) == 1
        models = [mdl.model for mdl in m.targets[0].models]
        assert models == ["af3"]
        runs = m.targets[0].models[0].runs
        assert len(runs) == 1
        assert runs[0].run_id == "200"
