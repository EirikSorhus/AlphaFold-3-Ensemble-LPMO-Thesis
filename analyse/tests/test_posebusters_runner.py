from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from lpmo_pipeline.qc.posebusters_runner import run_posebusters_single


def test_run_posebusters_single_auto_splits_combined_af3_export(tmp_path: Path) -> None:
    combined_pdb = tmp_path / "for_posebusters.pdb"
    combined_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY A   1      10.000  10.000  10.000  1.00  0.00           N",
                "ATOM      2  CA  GLY A   1      11.200  10.000  10.000  1.00  0.00           C",
                "HETATM    3  C1  NAG B   1      20.000  10.000  10.000  1.00  0.00           C",
                "HETATM    4  O4  NAG B   1      21.300  10.000  10.000  1.00  0.00           O",
                "HETATM    5  C1  NAG B   2      22.500  10.000  10.000  1.00  0.00           C",
                "CONECT    3    4    5",
                "CONECT    4    3",
                "CONECT    5    3",
                "END",
            ]
        )
        + "\n"
    )

    captured: dict[str, str | Path | None] = {}

    def _fake_run_pb_api(
        pdb_path: Path,
        protein_path: Path | None,
        reference_path: Path | None,
    ) -> dict[str, bool]:
        captured["ligand_path"] = pdb_path
        captured["protein_path"] = protein_path
        captured["reference_path"] = reference_path
        captured["ligand_text"] = pdb_path.read_text()
        captured["protein_text"] = protein_path.read_text() if protein_path is not None else None
        return {
            "all_atoms_connected": True,
            "internal_steric_clash": True,
            "inchi_convertible": True,
        }

    with patch("lpmo_pipeline.qc.posebusters_runner._run_pb_api", side_effect=_fake_run_pb_api):
        result = run_posebusters_single(combined_pdb, pose_id="pose_auto_split")

    assert result.passed is True
    assert captured["ligand_path"] != combined_pdb
    assert captured["protein_path"] is not None
    assert captured["reference_path"] is None

    ligand_text = str(captured["ligand_text"])
    protein_text = str(captured["protein_text"])

    assert " GLY A " not in ligand_text
    assert " NAG B " in ligand_text
    assert "CONECT    3    4    5" in ligand_text

    assert " GLY A " in protein_text
    assert " NAG B " not in protein_text
