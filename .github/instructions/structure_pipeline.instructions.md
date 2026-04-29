---
applyTo: "Masteroppgave/structure_pipeline/**"
---

# structure_pipeline – Repo-instruksjoner for Copilot og utviklere

> Ingen filendringer via kommandolinjeoperasjoner.

## Purpose
- `structure_pipeline/` er den manifestdrevne prediksjonspipelinen for AF3, Boltz-2 og RF3.
- Denne modulen er relativt stabil; foretrekk konservative feilrettinger og konfigurasjonsnare endringer fremfor redesign.
- `analyse/` skal konsumere artefakter herfra, ikke kopiere runnerlogikk.

## Read First
- `README.md` for faktisk CLI-arbeidsflyt, arkitektur og katalogkonvensjoner.
- `CODE_WALKTHROUGH.md` for praktisk orientering i kodebasen.
- `src/structure_pipeline/cli.py` for faktisk kommandooverflate.
- `src/structure_pipeline/config.py` for gjeldende konfigopplosning og path-regler.
- Runnerfiler under `src/structure_pipeline/runners/` nar du jobber naer modellspesifikk oppforsel.

## Active Conventions
- Pipelinen er manifestdrevet: `manifest/*.csv` og `DONE.ok` styrer resume og idempotens.
- AF3 bruker tostegs fan-in med separate MSA- og inferensjobber.
- `latest`-symlenker og `status.jsonl` brukes som operativ status, ikke bare logger.
- Relative input/output-stier kan losees via konfigfilen, men compute-node paths i path-seksjoner ma vaere gyldige for clusteret.
- Oligosakkaridspesifikasjon for AF3 horer hjemme i `config/oligo_definitions.yaml` og `runners/oligo.py`.

## Working Rules
- Ikke endre runner-semantikk uten a forsta hvordan den paverker manifestgenerering, jobbnavn, `DONE.ok`, `latest` og resume.
- Bevar tydelig skille mellom modellagnostisk orkestrering og modellspesifikk runnerkode.
- Unnga nye hardkodede cluster-paths i kildekoden; bruk konfigmodellene.
- Ved brukerendringer i CLI eller arbeidsflyt, oppdater relevant README-dokumentasjon.

## Validation
- For `python`/`pytest` i terminal: eksporter riktig env-bin i `PATH` en gang per terminalsession forst.
- Foretrekk smale tester i `tests/` og lette CLI-nare valideringer fremfor tunge jobbkjoringer.
- Bruk `--dry-run` eller tilsvarende trygg validering naer endringen berorer submit/resume-logikk.
- Ikke slett `DONE.ok` eller bruk tvungen rekjoring som del av vanlig validering uten eksplisitt behov.
- Små bash-kommandoer og svært små scripts kan kjores pa login node.
- Alle storre/tyngre oppgaver skal via SLURM, og agent skal ikke sende inn/styre SLURM-jobber uten eksplisitt brukertillatelse i chatten.

## Do
- Endre konfigmodeller, manifester og runnerkode i sma, sammenhengende steg.
- Hold FASTA-, ligand- og MSA-navnekonvensjoner stabile med mindre hele kjeden oppdateres.
- Legg til tester for nye kanttilfeller i parsing, config eller resumeoppforsel.

## Don't
- Ikke redesign hele pipelinen fra instruksjonsfiler eller spekulativ dokumentasjon alene.
- Ikke anta at lokal kjoring er lik SLURM-kjoring for tunge arbeidsflyter.
- Ikke introduser endringer som gjor `analyse/` avhengig av interne detaljer utover dokumenterte artefakter og CLI-kontrakter.
