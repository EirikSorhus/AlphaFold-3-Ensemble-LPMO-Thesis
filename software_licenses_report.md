# Programvare, pakker og lisenser brukt under `Masteroppgave_clean`

Dato for kartlegging: 2026-06-01.

Denne rapporten oppsummerer programvare, Python-pakker, eksterne programmer,
containere og dataressurser som brukes i prosjektmappene `analyse`,
`structure_pipeline`, `ligands` og `pipe_test`. Oversikten er basert på
`pyproject.toml`, `env.yml`-filene under `/cluster/work/projects/nn1003k/eirik/conda`,
pipelinekonfigurasjoner, shell-script og pakkemetadata hentet fra de
containeriserte Python-miljøene der dette lot seg gjøre.

Lisensfeltet er hentet fra pakkens egen Python-metadata når pakken er installert.
Der Conda- eller containerwrapperen ikke eksponerer full pakkemetadata, er
lisens markert som ikke verifisert fra lokalt miljø. For proprietære modeller,
vekter, databaser og eksterne webtjenester må lisens- og bruksvilkår håndteres
separat fra Python-pakkelisensene.

## Kort tekst til forskningsrapport

Analysene ble kjørt i isolerte, containeriserte Conda-miljøer på NRIS/HPC.
Strukturprediksjon ble orkestrert med en lokal `structure-pipeline` (MIT-lisens)
som sender SLURM-jobber til AlphaFold 3, Boltz-2 og RoseTTAFold3/Foundry via
Apptainer/Singularity-containere. Ligander ble representert som CCD/mmCIF,
MOL/SDF, PDB og SMILES, og konvertering/validering brukte blant annet RDKit og
gemmi. Nedstrøms struktur- og kontaktanalyse ble utført med en lokal
`lpmo-pipeline` (MIT-lisens) og tredjepartsverktøy som gemmi, MDAnalysis,
RDKit, ProLIF, PoseBusters, Privateer, HDBSCAN, scikit-learn, NumPy, SciPy og
pandas. Sekvensinnhenting og forbehandling i `pipe_test` brukte UniProt, CAZy,
NCBI/BLAST, InterPro, Biopython, requests og SignalP6, med domeneanalyser
planlagt mot dbCAN-HMM via HMMER. Full pakke- og lisensoversikt er gitt i
appendikset.

Merk: ProLIF rapporterer `0.0.0` i installert pakkemetadata i
`analyse_full_prolif_env`. Dette er ikke en brukbar versjonsverdi; riktig
ProLIF-versjon bør settes inn manuelt fra prosjektets installasjonslogg.

## 1. Analysepipeline (`Masteroppgave_clean/analyse`)

### Lokalt prosjekt

| Komponent | Versjon | Lisens | Rolle |
|---|---:|---|---|
| `lpmo-pipeline` | 0.1.0 | MIT | Lokal analysepipeline for QC, geometri, ProLIF-IFP, clustering, postprosessering og rapportering. |

### Miljø og programmer

Hovedmiljø: `/cluster/work/projects/nn1003k/eirik/conda/analyse_full_prolif_env`.

| Program/ressurs | Versjon/status | Lisens/status | Rolle |
|---|---:|---|---|
| Python | 3.11 | Python Software Foundation License | Kjøring av analysepipeline. |
| MAFFT | v7.526 | Ikke hentet fra lokalt metadatauttrekk | Sekvens-/familiejustering der brukt. |
| AmberTools-komponenter (`antechamber`, `sander`, `MMPBSA.py`, `pdb4amber`, m.fl.) | se pakkeoversikt | GPL/LGPL eller ikke oppgitt per komponent | Molekylforberedelse og kjemiverktøy i miljøet. Direkte `antechamber -h` svarte, men wrapperen manglet `amber.sh` i forventet sti. |
| Privateer container | `/cluster/projects/nn1003k/prog/privateer/privateer.sif` | Må verifiseres mot installert container/prosjektlisens | Glykan-/karbohydratvalidering. |
| PoseBusters container | `/cluster/projects/nn1003k/prog/posebusters/.../posebusters.sif` | PoseBusters Python-pakke: BSD License | Pose-QC og strukturell plausibilitet. |
| Apptainer/Singularity | systemverktøy | System-/distribusjonslisens, ikke hentet her | Kjøring av eksterne containere. |

### Direkte prosjektavhengigheter fra `analyse/pyproject.toml`

| Pakke/program | Minimum/krav | Lisens i installert metadata | Rolle |
|---|---:|---|---|
| gemmi | >=0.6.4 | MPL 2.0 | mmCIF/PDB-parsing og strukturhåndtering. |
| PoseBusters | >=0.3.0 | BSD License | Posevalidering. |
| Privateer | system/container | ikke verifisert lokalt | Karbohydrat-/glykanvalidering. |
| openbabel-wheel | >=3.1.1 | GPL | Molekylkonvertering/protonering. |
| reduce | AmberTools/MolProbity/system | ikke verifisert lokalt | Protonering/hydrogenplassering der brukt. |
| MDAnalysis | >=2.6.0 | LGPLv3+ | Struktur- og trajektorilignende analyse. |
| ProLIF | >=2.2.0 | Apache-2.0 i metadata, versjon rapporterer `0.0.0` | Protein-ligand interaction fingerprints. |
| hdbscan | >=0.8.33 | OSI Approved/BSD-lignende metadata | Tetthetsbasert clustering. |
| scipy | >=1.11.0 | BSD License | Numerikk/statistikk. |
| numpy | >=1.24.0 | BSD-3-Clause m.fl. | Numerikk. |
| pandas | >=2.0.0 | BSD License | Tabellbehandling. |
| scikit-learn | >=1.3.0 | BSD-3-Clause | Clustering, modeller og metrikker. |
| PyYAML | >=6.0 | MIT | Konfigurasjonsfiler. |
| jsonschema | >=4.20 | MIT | Validering av JSON-schema. |
| Jinja2 | >=3.1.2 | BSD License | Rapport-/templategenerering. |
| click | >=8.1.0 | BSD-3-Clause | CLI. |
| pytest, pytest-cov, ruff, mypy | dev | MIT/BSD/Apache, ruff/mypy ikke installert i miljøuttrekket | Testing og utvikling. |

### Full installert Python-pakkeoversikt for analyse

| Pakke | Versjon | Lisens fra metadata |
|---|---:|---|
| AmberUtils | 21.0 | GPL v2 or later |
| attrs | 26.1.0 | MIT |
| backports.zstd | 1.5.0 | PSF-2.0 |
| biopython | 1.87 | LicenseRef-Biopython-License-Agreement |
| Brotli | 1.2.0 | MIT |
| cached-property | 1.5.2 | BSD License |
| certifi | 2026.5.20 | Mozilla Public License 2.0 (MPL 2.0) |
| cffi | 2.0.0 | MIT |
| cftime | 1.6.5 | MIT |
| charset-normalizer | 3.4.7 | MIT |
| click | 8.4.0 | BSD-3-Clause |
| colorama | 0.4.6 | BSD License |
| contourpy | 1.3.3 | BSD License |
| coverage | 7.14.0 | Apache-2.0 |
| cycler | 0.12.1 | BSD License |
| dill | 0.4.1 | BSD License |
| edgembar | 3.0 | MIT |
| exceptiongroup | 1.3.1 | MIT License |
| filelock | 3.29.0 | MIT License |
| fonttools | 4.63.0 | MIT |
| freetype-py | 2.3.0 | BSD License |
| gemmi | 0.7.5 | Mozilla Public License 2.0 (MPL 2.0) |
| greenlet | 3.5.1 | MIT AND PSF-2.0 |
| GridDataFormats | 1.1.0 | LGPL-3.0-or-later |
| gsd | 5.0.1 | BSD-2-Clause |
| h2 | 4.3.0 | MIT License |
| h5py | 3.16.0 | BSD-3-Clause |
| hdbscan | 0.8.43 | OSI Approved |
| hpack | 4.1.0 | MIT License |
| hyperframe | 6.1.0 | MIT License |
| idna | 3.15 | BSD-3-Clause |
| iniconfig | 2.3.0 | MIT |
| Jinja2 | 3.1.6 | BSD License |
| joblib | 1.5.3 | BSD-3-Clause |
| jsonschema | 4.26.0 | MIT |
| jsonschema-specifications | 2025.9.1 | MIT |
| kiwisolver | 1.5.0 | BSD License |
| lpmo-pipeline | 0.1.0 | MIT |
| MarkupSafe | 3.0.3 | BSD-3-Clause |
| matplotlib | 3.10.9 | Python Software Foundation License |
| mda_xdrlib | 0.2.0 | Python Software Foundation License |
| MDAnalysis | 2.10.0 | GNU Lesser General Public License v3 or later (LGPLv3+) |
| mmcif-pdbx | 2.0.1 | CC0 1.0 Universal (CC0 1.0) Public Domain Dedication |
| MMPBSA.py | 16.0 | GPL v2 or later |
| mmtf-python | 1.1.3 | MIT License |
| mrcfile | 1.5.4 | BSD License |
| msgpack | 1.1.2 | Apache-2.0 |
| multiprocess | 0.70.19 | BSD License |
| munkres | 1.1.4 | Apache Software License |
| ndfes | 3.0 | MIT |
| netCDF4 | 1.7.2 | MIT License |
| networkx | 3.6.1 | BSD-3-Clause |
| numpy | 2.4.6 | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 |
| openbabel-wheel | 3.1.1.23 | GNU General Public License (GPL) |
| OpenMM | 8.5.1 | Python Software Foundation License (BSD-like) |
| packaging | 26.2 | Apache-2.0 OR BSD-2-Clause |
| packmol_memgen | 2025.1.29 | GNU General Public License v2 (GPLv2) |
| pandas | 3.0.3 | BSD License |
| ParmEd | 4.3.1 | LGPL |
| patsy | 1.0.2 | BSD License |
| pdb2pqr | 3.6.1 | BSD License |
| pdb4amber | 22.0 | Ikke oppgitt i pakkemetadata |
| pdbfixer | 1.12.0 | MIT License |
| pillow | 12.2.0 | MIT-CMU |
| pip | 26.1.1 | MIT |
| pluggy | 1.6.0 | MIT License |
| ply | 3.11 | BSD |
| Pmw | 2.1.1 | BSD License |
| posebusters | 0.6.5 | BSD License |
| prolif | 0.0.0 | Apache-2.0 |
| propka | 3.5.1 | GNU Lesser General Public License v2 (LGPLv2) |
| psutil | 7.2.2 | BSD-3-Clause |
| pycairo | 1.29.0 | LGPL-2.1-only OR MPL-1.1 |
| pycollada | 0.9.3 | BSD |
| pycparser | 3.0 | BSD-3-Clause |
| pyedr | 0.8.0 | GNU Lesser General Public License v2 or later (LGPLv2+) |
| Pygments | 2.20.0 | BSD-2-Clause |
| pykerberos | 1.2.4 | Apache Software License |
| pymol | 2.6.2 | Ikke oppgitt i pakkemetadata |
| pyMSMT | 22.0 | GPL v3 or later |
| pyparsing | 3.3.2 | MIT |
| PyQt5 | 5.15.11 | GPL v3 |
| PyQt5_sip | 12.17.0 | BSD-2-Clause |
| PySocks | 1.7.1 | BSD |
| pytest | 9.0.3 | MIT |
| pytest-cov | 7.1.0 | MIT License |
| python-dateutil | 2.9.0.post0 | BSD License; Apache Software License |
| pytng | 0.0.0 | BSD License |
| pytraj | 2.0.6 | GNU General Public License v3 (GPLv3) |
| PyYAML | 6.0.3 | MIT License |
| rdkit | 2025.3.6 | BSD-3-Clause |
| referencing | 0.37.0 | MIT |
| reportlab | 4.5.1 | BSD License |
| requests | 2.34.2 | Apache Software License |
| rlPyCairo | 0.4.0 | BSD License |
| rpds-py | 0.30.0 | MIT |
| sander | 22.0 | GPL v2 or later |
| scikit-learn | 1.8.0 | BSD-3-Clause |
| scikit-umfpack | 0.4.2 | BSD License |
| scipy | 1.17.1 | BSD License |
| seaborn | 0.13.2 | BSD License |
| setuptools | 82.0.1 | MIT |
| sip | 6.10.0 | BSD License |
| six | 1.17.0 | MIT License |
| SQLAlchemy | 2.0.49 | MIT |
| statsmodels | 0.14.6 | BSD License |
| threadpoolctl | 3.6.0 | BSD License |
| tidynamics | 1.1.2 | BSD License |
| toml | 0.10.2 | MIT License |
| tomli | 2.4.1 | MIT |
| tqdm | 4.67.3 | MPL-2.0 AND MIT |
| typing_extensions | 4.15.0 | PSF-2.0 |
| unicodedata2 | 17.0.1 | Apache License 2.0 |
| urllib3 | 2.7.0 | MIT |
| wheel | 0.47.0 | MIT |
| zstandard | 0.25.0 | BSD-3-Clause |

## 2. Structure pipeline (`Masteroppgave_clean/structure_pipeline`)

### Lokalt prosjekt

| Komponent | Versjon | Lisens | Rolle |
|---|---:|---|---|
| `structure-pipeline` | 0.1.0 | MIT | Lokal orkestrering av strukturprediksjon med AF3, Boltz-2 og RF3 på SLURM/HPC. |

### Miljø og eksterne modeller

Hovedmiljø: `/cluster/work/projects/nn1003k/eirik/conda/structure_pipeline_env`.

| Program/ressurs | Lokal sti/status | Lisens/status | Rolle |
|---|---|---|---|
| Python | 3.11 | Python Software Foundation License | Pipeline-CLI og jobbgenerering. |
| AlphaFold 3 | `/cluster/projects/nn1003k/prog/af3`, `af3_cpu_amd64.sif`, `af3_gpu_arm64.sif`, weights | Lisens/vilkår må verifiseres fra AlphaFold 3-installasjonen og tilgangsavtalen | Strukturprediksjon; AF3 genererer MSA med jackhmmer/nhmmer. |
| AlphaFold databaser | `/cluster/work/shared/alphafold_uncompressed.squashfs` | Databasevilkår må verifiseres separat | Sekvens-/template-databaser for AF3. |
| Boltz-2 | `/cluster/projects/nn1003k/prog/boltz/boltz2_alt.sif`, weights | Lisens/vilkår må verifiseres fra Boltz-2-installasjonen | Strukturprediksjon med preberegnet MSA. |
| RoseTTAFold3/Foundry | `/cluster/projects/nn1003k/prog/foundry/foundry.sif`, RF3 checkpoint | Lisens/vilkår må verifiseres fra Foundry/RF3-installasjonen | Strukturprediksjon med preberegnet MSA. |
| Apptainer/Singularity | systemverktøy | System-/distribusjonslisens, ikke hentet her | Containerkjøring. |
| SLURM | systemverktøy | System-/distribusjonslisens, ikke hentet her | Jobbkø og ressursstyring. |
| HMMER i AF3-container | `jackhmmer`, `nhmmer`, `hmmalign`, `hmmsearch`, `hmmbuild` inne i AF3-image | Ikke hentet fra containeren | AF3 MSA-steg. |

### Direkte prosjektavhengigheter fra `structure_pipeline/pyproject.toml`

| Pakke | Krav | Lisens i installert metadata | Rolle |
|---|---:|---|---|
| typer | >=0.9.0 | MIT | CLI. |
| pydantic | >=2.0 | MIT | Konfigurasjonsvalidering. |
| PyYAML | >=6.0 | MIT | YAML-konfigurasjon. |
| rich | >=13.0 | MIT | Terminaltabeller og logging. |
| pytest, pytest-cov | dev | MIT | Testing. |

### Full installert Python-pakkeoversikt for structure pipeline

| Pakke | Versjon | Lisens fra metadata |
|---|---:|---|
| annotated-doc | 0.0.4 | MIT |
| annotated-types | 0.7.0 | MIT License |
| click | 8.3.1 | BSD-3-Clause |
| coverage | 7.14.0 | Apache-2.0 |
| iniconfig | 2.3.0 | MIT |
| lpmo-pipeline | 0.1.0 | MIT |
| markdown-it-py | 4.0.0 | MIT License |
| mdurl | 0.1.2 | MIT License |
| packaging | 26.0 | Apache-2.0 OR BSD-2-Clause |
| pip | 26.0 | MIT |
| pluggy | 1.6.0 | MIT License |
| pydantic | 2.12.5 | MIT |
| pydantic_core | 2.41.5 | MIT |
| Pygments | 2.19.2 | BSD License |
| pytest | 9.0.3 | MIT |
| pytest-cov | 7.1.0 | MIT License |
| PyYAML | 6.0.3 | MIT License |
| rich | 14.3.2 | MIT License |
| setuptools | 80.10.2 | MIT |
| shellingham | 1.5.4 | ISC License (ISCL) |
| structure-pipeline | 0.1.0 | MIT |
| typer | 0.21.1 | MIT |
| typing-inspection | 0.4.2 | MIT |
| typing_extensions | 4.15.0 | PSF-2.0 |
| wheel | 0.46.3 | MIT |

## 3. Ligander (`Masteroppgave_clean/ligands`)

### Innhold og format

Ligandmappen inneholder oligosakkaridrepresentasjoner for cellulose, amylose og
NAG/chitin i flere formater:

| Format/filtype | Bruk | Lisens/status |
|---|---|---|
| `.smiles` | Kjemisk strengrepresentasjon | Ingen separat lisensmetadata funnet i filene. |
| `.mol`/`.sdf` | Molekylstrukturinput | Ingen separat lisensmetadata funnet i filene. |
| `.pdb` | 3D-koordinater/minimerte strukturer | Ingen separat lisensmetadata funnet i filene. |
| `.cif` | CCD/mmCIF-komponenter for AF3/Boltz/RF3 | Ingen separat lisensmetadata funnet i filene. |
| AF3 JSON-eksempler | Eksempelinput til AlphaFold 3 | Lokale prosjektfiler. |
| Boltz YAML-eksempler | Eksempelinput til Boltz-2 | Lokale prosjektfiler. |

### Ligandverktøy

Miljø for CCD-konvertering: `/cluster/work/projects/nn1003k/eirik/conda/boltz_ccd_env`.
Dette miljøets `env.yml` oppgir Python 3.11, pip, gemmi >=0.6 og RDKit >=2023.09.
Wrappermappen eksponerte ikke en kjørbar Python ved kartleggingen, så nøyaktig
installert versjon/lisensmetadata ble ikke hentet direkte fra miljøet.

| Verktøy/pakke | Krav/status | Lisens/status | Rolle |
|---|---:|---|---|
| Python | 3.11 i `env.yml` | Python Software Foundation License | Kjøring av ligandverktøy. |
| gemmi | >=0.6 i `env.yml`; analyse-miljø har 0.7.5 | MPL 2.0 | Parsing av CCD/mmCIF. |
| RDKit | >=2023.09 i `env.yml`; analyse-miljø har 2025.3.6 | BSD-3-Clause | Molekylobjekter, validering og pickle-konvertering. |
| pickle | Python standardbibliotek | Python Software Foundation License | Serialisering av RDKit Mol for Boltz-kompatibilitet. |

## 4. `pipe_test`

### Miljø og eksterne tjenester

Hovedmiljø for modul 2 og deler av `pipe_test`:
`/cluster/work/projects/nn1003k/eirik/conda/lpmo_pipe_env`.

`lpmo_pipe_env/env.yml` oppgir Python 3.11 og HMMER, men wrapperens `bin`-mappe
eksponerte bare Python/pip ved kartlegging. HMMER brukes likevel eksplisitt i
`run_module_3.sh` via `hmmscan` mot dbCAN-HMM-databaser.

| Program/tjeneste | Versjon/status | Lisens/status | Rolle |
|---|---:|---|---|
| Python | 3.11 | Python Software Foundation License | Sekvens- og metadataforbehandling. |
| HMMER/hmmscan | Oppgitt i `env.yml`; ikke synlig i wrapper-`bin` | HMMER-lisens ikke verifisert lokalt | Domeneannotering mot dbCAN-HMM. |
| dbCAN-HMM | `data/domains/dbcan/dbCAN-HMMdb-V14.hmm`, `dbCAN_sub.hmm` | Databasevilkår må verifiseres separat | Karbohydrataktive enzymdomener. |
| SignalP6 | brukt via `$SIGNALP6_PATH` i `run_signalpeptide.sh` | SignalP6-lisens/vilkår må verifiseres separat | Signalpeptidprediksjon og trimming. |
| UniProt REST API | webtjeneste | UniProt-bruksvilkår/databaselisens må følges | Sekvens- og metadatahenting. |
| CAZy | webtjeneste/datafiler | CAZy-bruksvilkår må følges | Familieklassifisering og proteinlister. |
| NCBI Entrez/eUtils | webtjeneste | NCBI-bruksvilkår må følges | Fallback-henting av sekvenser/metadata. |
| NCBI BLAST | webtjeneste | NCBI BLAST-bruksvilkår må følges | Sekvenssøk mot Swiss-Prot ved ukjente ID-er. |
| InterPro API | webtjeneste/cache | InterPro-bruksvilkår må følges | Domene- og SignalP-relatert metadatafallback. |

### Full installert Python-pakkeoversikt for `lpmo_pipe_env`

| Pakke | Versjon | Lisens fra metadata |
|---|---:|---|
| appdirs | 1.4.4 | MIT License |
| attrs | 25.4.0 | MIT |
| beautifulsoup4 | 4.14.3 | MIT License |
| biopython | 1.86 | Freely Distributable |
| bioservices | 1.12.1 | BSD License |
| cattrs | 25.3.0 | MIT License |
| certifi | 2025.11.12 | Mozilla Public License 2.0 (MPL 2.0) |
| charset-normalizer | 3.4.4 | MIT |
| click | 8.3.1 | BSD-3-Clause |
| colorama | 0.4.6 | BSD License |
| colorlog | 6.10.1 | MIT License |
| contourpy | 1.3.3 | BSD License |
| coverage | 7.14.0 | Apache-2.0 |
| cycler | 0.12.1 | BSD License |
| easydev | 0.13.3 | BSD License |
| fonttools | 4.61.0 | MIT |
| gevent | 25.9.1 | MIT |
| greenlet | 3.3.0 | MIT AND Python-2.0 |
| grequests | 0.7.0 | BSD License |
| idna | 3.11 | BSD-3-Clause |
| iniconfig | 2.3.0 | MIT |
| kiwisolver | 1.4.9 | BSD License |
| line_profiler | 4.2.0 | BSD License |
| lpmo-pipeline | 0.1.0 | MIT |
| lxml | 5.4.0 | BSD License |
| markdown-it-py | 4.0.0 | MIT License |
| matplotlib | 3.10.8 | Python Software Foundation License |
| mdurl | 0.1.2 | MIT License |
| numpy | 2.3.5 | BSD License |
| packaging | 25.0 | Apache Software License; BSD License |
| pandas | 2.3.3 | BSD License |
| pexpect | 4.9.0 | ISC License (ISCL) |
| pillow | 12.0.0 | MIT-CMU |
| pip | 25.3 | MIT |
| platformdirs | 4.5.1 | MIT License |
| pluggy | 1.6.0 | MIT License |
| ptyprocess | 0.7.0 | ISC License (ISCL) |
| Pygments | 2.19.2 | BSD License |
| pyparsing | 3.2.5 | MIT |
| pytest | 9.0.3 | MIT |
| pytest-cov | 7.1.0 | MIT License |
| python-dateutil | 2.9.0.post0 | BSD License; Apache Software License |
| pytz | 2025.2 | MIT License |
| PyYAML | 6.0.3 | MIT License |
| requests | 2.32.5 | Apache Software License |
| requests-cache | 1.2.1 | BSD License |
| rich | 14.2.0 | MIT License |
| rich-click | 1.9.4 | MIT License |
| setuptools | 80.9.0 | MIT |
| six | 1.17.0 | MIT License |
| soupsieve | 2.8 | MIT License |
| suds-community | 1.2.0 | GNU Library or Lesser General Public License (LGPL) |
| tqdm | 4.67.1 | MIT License; Mozilla Public License 2.0 (MPL 2.0) |
| typing_extensions | 4.15.0 | PSF-2.0 |
| tzdata | 2025.2 | Apache Software License |
| url-normalize | 2.2.1 | MIT |
| urllib3 | 2.6.2 | MIT |
| wheel | 0.45.1 | MIT License |
| wrapt | 1.17.3 | BSD License |
| xmltodict | 0.14.2 | MIT License |
| zope.event | 6.1 | ZPL-2.1 |
| zope.interface | 8.1.1 | Zope Public License |

### SignalP6-miljø

Miljø: `/cluster/work/projects/nn1003k/eirik/conda/signalp6_env`.

Wrappermappen hadde ingen synlige kjørbare filer under kartleggingen.
`env.yml` oppgir:

| Pakke | Krav | Lisens/status |
|---|---:|---|
| Python | 3.10 | Python Software Foundation License |
| pip | ikke låst | MIT |
| PyTorch | ikke låst | BSD-style/PyTorch-lisens må verifiseres mot installasjon |
| NumPy | <2 | BSD |

`run_signalpeptide.sh` bruker i tillegg ekstern `signalp6`-binær fra
`$SIGNALP6_PATH` med standardsti `/cluster/home/eisorhus/.local/bin/signalp6`.
SignalP6-versjon og lisens må dokumenteres fra den installasjonen.

## 5. Andre miljøer under `/conda`

Følgende miljøer finnes i prosjektroten, men er ikke hovedmiljøene for den
rapporterte kjøringen eller har delvis overlapp:

| Miljø | Innhold fra `env.yml` | Status |
|---|---|---|
| `analyse_env` | Ligner analysemiljøet: Python 3.11, NumPy, AmberTools, gemmi, hdbscan, MDAnalysis, RDKit, scikit-learn, scipy, m.fl. | Eldre/alternativt analysemiljø. Ikke brukt som primærsti i `runtime_paths.yaml`. |
| `analyse_new_prolif_env` | Alternativt analyse-/ProLIF-miljø | Ikke kartlagt i detalj; hovedpipeline peker til `analyse_full_prolif_env`. |
| `boltz_ccd_env` | Python 3.11, pip, gemmi >=0.6, RDKit >=2023.09 | Brukt/relevant for ligandkonvertering. |
| `lpmo_pipe_env` | Python 3.11, HMMER | Brukt/relevant for `pipe_test`. |
| `signalp6_env` | Python 3.10, pip, PyTorch, NumPy<2 | Relevant for SignalP6-kjøring. |
| `structure_pipeline_env` | Python 3.11, PyYAML, Typer, Pydantic, Rich, pytest | Brukt for `structure_pipeline`. |

## 6. Lisensmessige forbehold

1. Python-pakkelisenser er hentet fra installert pakkemetadata. Dette er egnet
   for rapportering, men ikke en juridisk fullstendig lisensrevisjon.
2. Containere, modellvekter og databaser for AlphaFold 3, Boltz-2,
   RoseTTAFold3/Foundry, Privateer, SignalP6, dbCAN-HMM og AlphaFold-databaser
   må dokumenteres med egne lisens- og tilgangsvilkår fra installasjonsstedet.
3. Eksterne webtjenester som UniProt, CAZy, NCBI BLAST/Entrez og InterPro er
   datakilder/tjenester, ikke bare programvarepakker. Bruksvilkår og sitering
   bør oppgis separat i metodedel/referanseliste.
4. Ligandfilene i `Masteroppgave_clean/ligands` inneholder ikke innebygd
   lisensmetadata. Hvis de er generert fra eksterne databaser eller verktøy,
   bør opprinnelig kilde, verktøy og eventuelle vilkår dokumenteres.
5. ProLIF-versjonen må korrigeres manuelt, siden installert metadata rapporterer
   `0.0.0` selv om faktisk brukt versjon er kjent utenfor metadatafeltet.
