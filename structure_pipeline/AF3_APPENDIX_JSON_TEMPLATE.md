# Appendix template: AlphaFold 3 JSON input

The AlphaFold 3 input files were generated as JSON documents with the `alphafold3`
dialect. The example below is a shortened version of the CEL4 input file for
`A0A0A1ED04_CEL4`. Long text fields, such as the MSA strings and mmCIF template
coordinates, are abbreviated with placeholders to keep the appendix readable.

```json
{
  "dialect": "alphafold3",
  "version": 4,
  "name": "A0A0A1ED04_CEL4",
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "HGWIEESRAGLCMTGQNTGCGAVQYEPWSVEGRGDFPEIGVPDGEITGGGKYAPLYEQTATRWTKVNMTGGPYTFHWKMVANHSTNRWDYYITKPGWNPNEPIGRDDIELFCRYEDNGAIPPMDVLNDCYIPNDREGYHVIIGVWDIFDTVNAFYQAIDV",
        "modifications": [],
        "unpairedMsa": ">query\nHGWIEESRAGLCMTGQNTGCGAVQYEPWSVEGRGDFPEIGVPDGEITGGGKYAPLYEQTATRWTKVNMTGGPYTFHWKMVANHSTNRWDYYITKPGWNPNEPIGRDDIELFCRYEDNGAIPPMDVLNDCYIPNDREGYHVIIGVWDIFDTVNAFYQAIDV\n>sequence_hit_1\n<aligned_sequence_1>\n>sequence_hit_2\n<aligned_sequence_2>\n...",
        "pairedMsa": ">query\nHGWIEESRAGLCMTGQNTGCGAVQYEPWSVEGRGDFPEIGVPDGEITGGGKYAPLYEQTATRWTKVNMTGGPYTFHWKMVANHSTNRWDYYITKPGWNPNEPIGRDDIELFCRYEDNGAIPPMDVLNDCYIPNDREGYHVIIGVWDIFDTVNAFYQAIDV\n>paired_hit_1\n<aligned_sequence_1>\n>paired_hit_2\n<aligned_sequence_2>\n...",
        "templates": [
          {
            "mmcif": "data_template_1\n#\n<mmCIF template coordinates and metadata>\n#",
            "queryIndices": [0, 1, 2, 3, 4, 5],
            "templateIndices": [0, 1, 2, 3, 4, 5]
          }
        ]
      }
    },
    {
      "ligand": {
        "id": "B",
        "ccdCodes": ["CU"]
      }
    },
    {
      "ligand": {
        "id": "C",
        "ccdCodes": ["BGC", "BGC", "BGC", "BGC"]
      }
    }
  ],
  "modelSeeds": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
  "bondedAtomPairs": [
    [["C", 1, "C1"], ["C", 2, "O4"]],
    [["C", 2, "C1"], ["C", 3, "O4"]],
    [["C", 3, "C1"], ["C", 4, "O4"]]
  ]
}
```

In the full generated files, `unpairedMsa` and `pairedMsa` contain complete FASTA
formatted multiple sequence alignments. The `templates` list contains one object
per structural template; in this CEL4 example the full file contains four
template objects. Each template object stores the template as an mmCIF string and
the matching zero-based residue indices in the query and template sequences. For
this protein, each full index list contains 160 aligned residues. The two ligand
entries define the catalytic copper ion (`CU`) and the CEL4 substrate,
represented as four linked beta-D-glucose CCD residues (`BGC`). The
`bondedAtomPairs` entries define the glycosidic bonds between consecutive glucose
units in ligand chain `C`.
