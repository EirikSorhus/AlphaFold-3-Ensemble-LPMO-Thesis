# Global Copilot Instructions

## Scope
- Gjelder hele repoet.
- Ikke finn pa verktoy, kommandoer, paths eller konfigurasjoner som ikke er dokumentert.
- Hvis informasjon mangler: bruk `TODO` eller `PLACEHOLDER`.
- Ingen endringer, slettinger eller laging av nye filer skal skje ved kommandolinje operasjoner.

## Project Context
- Repoet utvikler en HPC-basert pipeline for strukturprediksjon og analyse av LPMO-enzymer.
- Overordnet flyt: `metadata -> sequences -> ligand pairing -> structure prediction -> parsing -> analysis`.
- Strukturprediksjon bruker `AF3`, `RF3` og `Boltz-2`.

## HPC Rules (Critical)
- All kode som ikke er ekstremt lett skal kjores via SLURM.
- Anta at compute-jobs kjores med `sbatch` eller tilsvarende scheduler.
- Ikke implementer lokal kjaring for tunge oppgaver.
- Design pipelinekode for batch-jobber, ikke interaktiv kjaring.
- Gjor scripts kompatible med HPC job arrays nar relevant.

## Python and Environments
- System-Python pa HPC-clusteret er for gammel; ikke anta system-Python.
- All Python-kode skal kjores i aktivert conda-environment.
- Environmenter ligger i `eirik/conda/`.
- Anta at riktig environment er aktivert forst.
- Anta aldri environment uten dokumentasjon.

## Required Submodule Instructions
- Les `.github/instructions/<submodule>.md` for aktuell mappe for du foreslar endringer.
- Hver slik fil skal beskrive: formal, conda-environment, environment-aktivering, HPC-kjoring.
- Hvis fil mangler eller er uklar: bruk `TODO` eller `PLACEHOLDER` i stedet for antakelser.

## Repository Priorities
- `ligands/`: ligandfiler og filkonvertering.
- `structure_pipeline/`: i stor grad ferdig; kun sma, konservative endringer.
- `analyse/`: hovedsakelig pseudokode; mangler implementasjon og logikkgjennomgang.
- `pipe_test/`: fokuser kun pa `pipe_test/scripts/module_2_new/` og `pipe_test/data/`.
- Andre `module_2*`-varianter i `pipe_test/scripts/` er gamle/ufullstendige og skal normalt ignoreres.

## Frozen Code
- `pipe_test/scripts/module_2_new/` er stabil, ferdig og behandles som frozen code.
- Ikke refaktorer eller endre denne modulen uten eksplisitt brukerinstruksjon.
- Hvis endringer i frozen code foreslas: merk dem tydelig som eksplisitt avvik.
- Shell scripts som kjorer `module_2_new` skal ikke brukes eller endres uten eksplisitt foresporsel.

## Change Policy
- Foretrekk sma, sikre og isolerte endringer.
- Ikke endre mange moduler samtidig uten klar grunn.
- Ikke bryt eksisterende IO-kontrakter.
- Ikke fjern hardkodede paths uten dokumentert arsak (kan vaere nodvendig pa HPC).
- I `analyse/`: vaer tydelig pa hva som er ny implementasjon vs eksisterende pseudokode.

## README Rule
- Oppdater `README.md` hver gang kode endres.
- Dokumenter nye scripts, workflow-endringer, nye steg og nye avhengigheter.
- Foresla alltid README-endringer nar kode endres.

## Pipeline Design Principles
- Hold pipeline-steg tydelig separert.
- Gjor IO eksplisitt med tydelige input/output paths.
- Unnga unodvendig kopiering av store filer.

## Logging Requirements
- Alle pipeline-steg bor logge inputfiler.
- Alle pipeline-steg bor logge output paths.
- Alle pipeline-steg bor logge parametere.
- Alle pipeline-steg bor logge timestamps.
- Alle pipeline-steg bor logge verktoy som brukes (`AF3`, `RF3`, `Boltz-2`).

## Error Handling
- Fail fast ved manglende input.
- Feilmeldinger skal forklare hva som feilet.
- Feilmeldinger skal forklare hvor det feilet.
- Feilmeldinger skal forklare hvordan det kan fikses.

## Long-Term Configuration Direction
- Layered configuration management er et langsiktig maal.
- Malet er multi-tool HPC pipelines med lagdelte config-filer.
- Nar relevant kan struktur foreslas (for eksempel global + tool + run config).
- Ikke implementer full konfigurasjonsplattform uten eksplisitt instruksjon.