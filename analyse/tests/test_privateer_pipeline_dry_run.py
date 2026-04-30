from __future__ import annotations

from pathlib import Path

from lpmo_pipeline.qc.privateer_runner import prepare_privateer_input, run_privateer


def test_privateer_pipeline_dry_run_covers_prep_and_command_build(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from lpmo_pipeline.qc import privateer_runner as mod

    fake_sif = tmp_path / "privateer.sif"
    fake_sif.write_text("fake")
    monkeypatch.setattr(mod, "PRIVATEER_SIF_CANDIDATES", (fake_sif,))

    source = tmp_path / "input.cif"
    source.write_text(
        "\n".join(
            [
                "data_test",
                "#",
                "loop_",
                "_atom_site.group_PDB",
                "_atom_site.id",
                "_atom_site.type_symbol",
                "_atom_site.label_atom_id",
                "_atom_site.label_alt_id",
                "_atom_site.label_comp_id",
                "_atom_site.label_asym_id",
                "_atom_site.label_entity_id",
                "_atom_site.label_seq_id",
                "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.occupancy",
                "_atom_site.B_iso_or_equiv",
                "_atom_site.auth_seq_id",
                "_atom_site.auth_asym_id",
                "_atom_site.pdbx_PDB_model_num",
                "HETATM 1 C C1 . NAG C 1 1 ? 0.000 0.000 0.000 1.00 20.00 1 C 1",
                "HETATM 2 O O5 . NAG C 1 1 ? 1.200 0.000 0.000 1.00 20.00 1 C 1",
                "#",
                "",
            ]
        )
    )

    prepared = prepare_privateer_input(source, tmp_path / "privateer_input.cif")
    result = run_privateer(
        prepared,
        pose_id="pose dry run whole slice",
        output_dir=tmp_path / "privateer_work",
        dry_run=True,
    )

    assert prepared.exists()
    assert result.dry_run is True
    assert result.command[:3] == ["apptainer", "run", "--cleanenv"]
    assert "-pdbin" in result.command
    assert str(prepared.resolve()) in result.command
    assert result.validation_data_path.endswith("validation_data-privateer")