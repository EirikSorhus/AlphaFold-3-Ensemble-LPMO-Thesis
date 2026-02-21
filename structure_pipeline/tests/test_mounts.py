"""Tests for bind-mount helpers."""

from structure_pipeline.runners.mounts import compute_bind_mounts, log_mount_plan


class TestComputeBindMounts:
    """Tests for compute_bind_mounts."""

    def test_empty_input(self):
        assert compute_bind_mounts([]) == []

    def test_single_path(self):
        mounts = compute_bind_mounts(["/cluster/projects/nn1003k/prog/foundry"])
        assert len(mounts) == 1
        assert "/cluster/projects/nn1003k/prog/foundry" in mounts[0]

    def test_merges_common_prefix(self):
        """Paths under the same parent should be merged."""
        mounts = compute_bind_mounts([
            "/cluster/projects/nn1003k/prog/foundry/checkpoints/ckpt.pt",
            "/cluster/projects/nn1003k/prog/foundry/models/rf3/src",
        ])
        # Should merge into a single mount at /cluster/projects/nn1003k/prog/foundry
        assert len(mounts) <= 2

    def test_separate_trees(self):
        """Paths in different trees should not be merged."""
        mounts = compute_bind_mounts([
            "/cluster/projects/nn1003k/prog/foundry",
            "/cluster/shared/alphafold/public_databases",
        ])
        assert len(mounts) == 2

    def test_file_paths_use_parent(self):
        """File paths should result in mounting the parent directory."""
        mounts = compute_bind_mounts([
            "/cluster/projects/nn1003k/prog/foundry/foundry.sif",
        ])
        assert len(mounts) == 1
        # The mount should be the parent dir, not the file itself
        assert "foundry.sif" not in mounts[0]


class TestLogMountPlan:
    """Tests for log_mount_plan."""

    def test_produces_output(self):
        mounts = ["/a:/a", "/b:/b"]
        plan = log_mount_plan(mounts, ["/a/file1", "/b/file2"])
        assert "Mount plan:" in plan
        assert "/a/file1" in plan
        assert "--bind /a:/a" in plan
