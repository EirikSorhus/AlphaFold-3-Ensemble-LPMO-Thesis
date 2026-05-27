from __future__ import annotations

from pathlib import Path

from lpmo_pipeline.analysis.residue_region_annotation import (
    build_region_annotations_for_ifp_results,
    load_protein_region_definitions,
)
from lpmo_pipeline.analysis.prolif_ifp import IFPResult


def test_metadata_core_range_is_mature_relative(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "UniProt_ID\tSignal_End\tLPMO_Core_Start\tLPMO_Core_End\n"
        "P1\t20\t21\t120\n"
    )

    definitions = load_protein_region_definitions(metadata_path)

    assert definitions["P1"].core_start == 1
    assert definitions["P1"].core_end == 100


def test_region_annotations_mark_core_and_non_core_for_full_length(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "UniProt_ID\tSignal_End\tLPMO_Core_Start\tLPMO_Core_End\n"
        "P1\t20\t21\t120\n"
    )
    definitions = load_protein_region_definitions(metadata_path)
    result = IFPResult(
        pose_id="pose-1",
        feature_names=[
            "NAG1.B|ASN10.A|ImplicitHBDonor",
            "NAG1.B|TYR130.A|ImplicitHBAcceptor",
        ],
        flat_bitvector=[1, 1],
    )

    annotations = build_region_annotations_for_ifp_results(
        [result],
        protein_id="P1",
        construct_type="full_length",
        region_definitions=definitions,
    )

    assert annotations["pose-1"][("A", 10, "ASN")].is_core_region is True
    assert annotations["pose-1"][("A", 10, "ASN")].is_non_core_region is False
    assert annotations["pose-1"][("A", 130, "TYR")].is_core_region is False
    assert annotations["pose-1"][("A", 130, "TYR")].is_non_core_region is True
