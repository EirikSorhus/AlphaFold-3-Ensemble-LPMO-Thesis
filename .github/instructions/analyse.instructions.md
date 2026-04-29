---
applyTo: "Masteroppgave/analyse/**"
description: "Use for analyse/ work: AF3-only LPMO analysis pipeline, governing document priority, hard invariants, and focused validation commands."
---

# analyse – Repo-instruksjoner for Copilot og utviklere

> Status: aktiv Python-pakke med paagaaende implementasjon. Ikke anta at eldre pseudokode eller arkiverte planer fortsatt gjelder.
> Ingen endringer, slettinger eller laging av nye filer skal skje ved kommandolinje operasjoner.

## Purpose
- `analyse/` er analyse- og rapportlaget for LPMO-prosjektet.
- Aktivt scope er AF3-only analyse av predikerte strukturer; RF3 og Boltz-2 er ekskludert fra analysearbeidsflyten.
- Koden skal konsumere prediksjonsartefakter fra `structure_pipeline/` og ikke duplisere modellrunner-logikk.

## Read First
- `README.md` for aktiv status, quick start og dokumentasjonskart.
- `AF3_LPMO_pipeline_detailed_plan.md` som primar styringskilde for analysedesign.
- `MASTERPLAN.md` og `IMPLEMENTATION_PLAYBOOK.md` for implementasjonsrekkefolge, artefaktkrav og stopppunkter.
- `OPEN_QUESTIONS.md` for avklarte og uavklarte beslutninger som ikke skal gjettes.
- `copilot.md` for prosjektspesifikke AI-regler og ansvarsskille mot `structure_pipeline/`.

## Working Rules
- Bevar skillet mellom prediksjon og analyse: nye behov for AF3-kjoring eller parserlogikk i prediksjonslaget horer hjemme i `structure_pipeline/`.
- Implementer en modul eller ett naerliggende artefaktsteg av gangen.
- Hold endringer kompatible med eksisterende artefaktkontrakter og skjemaer.
- Hvis du ma velge mellom motstridende dokumenter, folg prioriteten oppgitt i `README.md` og `MASTERPLAN.md`.

## Active Conventions
- Kanonisk term er `protein` og `protein_id` i ny kode.
- Hovedanalyse er cluster-primary; ikke kollaps cluster-rader til protein som primar analyseenhet.
- Hovedinput er AF3-artefakter i mmCIF-format.
- Kjedeoppsett i normaliserte strukturer er protein `A`, glykaner `B..D`, metall `E`.
- Pre-QC aktivt sete-proximity gate korer foran PoseBusters og Privateer.
- IFP-clustering bruker IFP-features; geometri legges pa etter clustering.
- Hardcoded site-spesifikke paths skal unngas i kildekode.

## Environment And Validation
- Bruk dokumentert analysemiljo under `/cluster/work/projects/nn1003k/eirik/conda/` for Python-kjoring.
- For `python`/`pytest` i terminal: eksporter riktig env-bin i `PATH` en gang per terminalsession forst.
- Ikke installer pakker som del av vanlig kodearbeid.
- Foretrekk smale tester eller kontraktvalidering fremfor tunge end-to-end-kjoringer.
- Hvis ekstern programvare som `privateer` eller `reduce` mangler, rapporter det i stedet for a improvisere lokale workarounds.
- Små bash-kommandoer og svært små scripts kan kjores pa login node; storre oppgaver skal via SLURM.
- Agent skal ikke sende inn/styre SLURM-jobber uten eksplisitt tillatelse fra bruker i chatten.

## Do
- Oppdater tester i `tests/` naer endret oppforsel.
- Valider skjema- eller artefaktendringer mot `schemas/` og eksisterende testkontrakter.
- Marker ekte uavklarte designvalg i kode eller dokumentasjon i stedet for a gjette.
- Lenke til styrende dokumentasjon i svar og nye customization-filer fremfor a kopiere store utdrag.

## Don't
- Ikke reimplementer AF3-, RF3- eller Boltz-runnerlogikk i `analyse/`.
- Ikke behandle PLACER som aktivt, obligatorisk steg uten eksplisitt dokumentert grunnlag i dagens planfiler.
- Ikke bruk `pipe_test/` som autoritativ kilde for analysekontrakter.
- Ikke endre låste analyseparametre eller operative terskler uten at styrende dokumentasjon eller brukerbeskjed sier det.
- Hold thresholds og tunables i config/schema, ikke som skjulte literals i kode.
