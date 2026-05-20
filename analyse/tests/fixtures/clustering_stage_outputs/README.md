Stage 16 clustering-output fixture for downstream cluster-dependent tests.

These files are small, deterministic Stage 16 surfaces:
- `cluster_signatures.json`
- `cluster_residue_signature.tsv`
- `cluster_ifp_signature.tsv`

Use this fixture for tests that require non-empty clusters without rerunning
the full real-data clustering pipeline. It is currently the canonical local
fixture for Stage 16b residue-importance tests because the short real-data
smoke usually contains too few analyzed poses to produce retained clusters.
