#!/usr/bin/env bash
# submit MSA, then inference once MSA is done
sbatch --dependency=afterok:$(sbatch --parsable af3_msa.sh) af3_inf.sh
