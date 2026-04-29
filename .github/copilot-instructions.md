# Masteroppgave Copilot Instructions

## Scope
- Gjelder hele `Masteroppgave/`.
- Ikke finn pa verktoy, paths, kommandoer eller miljovalg som ikke er dokumentert i repoet.
- Hvis informasjon mangler, bruk `TODO` eller `PLACEHOLDER` i stedet for antakelser.
- Ikke opprett, endre eller slett filer via shell-kommandoer; bruk editor-baserte filendringer.

## First Routing Rule
- Les aktuell filinstruksjon i `.github/instructions/` for delomradet du jobber i for du foreslar eller gjor endringer.
- Bruk disse som primarkilde for mappe-spesifikke regler:
  - `analyse/` -> `.github/instructions/analyse.instructions.md`
  - `structure_pipeline/` -> `.github/instructions/structure_pipeline.instructions.md`
  - `ligands/` -> `.github/instructions/ligands.instructions.md`
  - `pipe_test/scripts/module_2_new/` -> `.github/instructions/module_2_new.instructions.md`

## Project Map
- Repoet inneholder en HPC-orientert arbeidsflyt for LPMO-prosjekter.
- `structure_pipeline/` er produksjonsnart manifestdrevet prediksjonspipeline for `AF3`, `Boltz-2` og `RF3`.
- `analyse/` er analysepipen og er na AF3-only; den skal konsumere prediksjonsartefakter og ikke reimplementere runner-logikk.
- `ligands/` inneholder ligandfiler og konverteringshjelpere, ikke hoved-pipelineorkestrering.
- `pipe_test/` er delvis eksperimentell; kun `pipe_test/scripts/module_2_new/` behandles som stabilt arbeidsomrade.

## HPC And Environment Rules
- Anta HPC-kjoring med SLURM for alt som er mer enn lett lokal validering.
- Ikke design tunge workflows rundt interaktiv lokal kjoring.
- System-Python pa clusteret skal ikke antas brukbar; bruk dokumentert conda-environment under `/cluster/work/projects/nn1003k/eirik/conda/`.
- Hvis riktig environment ikke er dokumentert for delprosjektet, stopp og marker mangelen i stedet for a gjette.

## Terminal Session Rules
- For `python` eller `pytest`: eksporter riktig environment-bin i `PATH` for aktiv terminalsession forst (en gang per terminal), for eksempel `export PATH="/cluster/work/projects/nn1003k/eirik/conda/<env>/bin:$PATH"`.
- Små bash-kommandoer og svært små scripts kan kjores direkte pa login node.
- Alle storre/tyngre oppgaver skal kjores via SLURM.
- Kun bruker skal sende inn eller styre SLURM-jobber, med mindre brukeren gir eksplisitt tillatelse i chatten.

## Change Policy
- Foretrekk sma, isolerte endringer som bevarer eksisterende IO-kontrakter.
- Ikke endre hardkodede HPC-stier uten dokumentert grunnlag i repoet.
- Bevar separasjonen mellom prediksjon (`structure_pipeline/`) og analyse (`analyse/`).
- Naerliggende README eller annen styrende dokumentasjon skal oppdateres nar kodeendringen faktisk endrer bruk, input/output eller operative antakelser.

## Frozen Areas
- `pipe_test/scripts/module_2_new/` behandles som frozen code med mindre brukeren eksplisitt ber om endring der.
- Andre `module_2*`-varianter under `pipe_test/scripts/` er eldre og skal normalt ignoreres.

## Design And Reliability
- Hold pipeline-steg tydelig separert med eksplisitte input/output paths.
- Logg sentrale parametere, inputkilder, output paths og gate-utfall for pipeline-steg.
- Feilmeldinger skal vaere handlingsrettede: hva som feilet, hvor det feilet, og hva som ma sjekkes.

## Authoritative Documentation
- For `analyse/`, bruk README + planfilene som kilde for operative beslutninger; ikke stol pa eldre pseudokode alene.
- For `structure_pipeline/`, bruk `README.md`, `CODE_WALKTHROUGH.md`, CLI-en og konfigmodellene som autoritative beskrivelser av faktisk oppforsel.