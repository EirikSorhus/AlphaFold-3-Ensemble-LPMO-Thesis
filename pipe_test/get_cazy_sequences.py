import argparse
import requests
from datetime import datetime
import uuid
import json

# Parse input arguments for families
parser = argparse.ArgumentParser(description="Hent sekvenser basert på CAZy-familier (UniProt API)")
parser.add_argument("-f", "--families", nargs="+", required=True, help="Liste over CAZy-familier som skal hentes")
args = parser.parse_args()

families = args.families

# Opprett unikt run ID, tidsstempel og datostreng for filnavn
now = datetime.now()
timestamp_str = now.isoformat()
date_str = now.strftime("%Y%m%d")
run_id = str(uuid.uuid4())

# Forbered datastrukturer for å samle resultater
seen_ids = set()
combined_data = {}
family_counts = {}

# Åpne kombinert FASTA-fil på forhånd (slik at vi kan skrive fortløpende)
combined_fasta_filename = f"all_sequences_{date_str}.fasta"
combined_fasta_file = open(combined_fasta_filename, "w")

# Iterer over hver etterspurte familie og hent sekvenser fra UniProt
for fam in families:
    # Bygg spørrings-URL for UniProt REST API (bruk stream-endepunktet for å få alle resultater)
    query = f"xref:CAZY-{fam}"
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {"query": query, "format": "fasta"}
    print(f"[INFO] Henter sekvenser for familie {fam} ...")
    try:
        response = requests.get(url, params=params, stream=True)
        response.raise_for_status()
    except Exception as e:
        print(f"[ERROR] Feil ved henting av {fam}: {e}")
        family_counts[fam] = 0
        continue

    # Åpne egen FASTA-fil for denne familien
    family_fasta_filename = f"{fam}_sequences_{date_str}.fasta"
    fam_file = open(family_fasta_filename, "w")

    count = 0  # teller sekvenser for denne familien
    current_seq_id = None
    current_seq_name = None
    current_seq_org = None
    seq_buffer = []

    # Stream gjennom linjene i FASTA-responsen
    for line in response.iter_lines(decode_unicode=True):
        if not line:
            continue  # hopp over tomme linjer
        if line.startswith(">"):
            # Hvis vi allerede har en pågående sekvens, avslutt den først
            if current_seq_id is not None:
                seq_sequence = "".join(seq_buffer)
                seq_length = len(seq_sequence)
                # Skriv sekvensen til familie-FASTA
                fam_file.write(f">{current_seq_id} {current_seq_name} [{current_seq_org}]\n")
                for i in range(0, len(seq_sequence), 60):
                    fam_file.write(seq_sequence[i:i+60] + "\n")
                # Legg til i kombinert FASTA og metadata hvis ny sekvens
                if current_seq_id not in seen_ids:
                    combined_fasta_file.write(f">{current_seq_id} {current_seq_name} [{current_seq_org}]\n")
                    for i in range(0, len(seq_sequence), 60):
                        combined_fasta_file.write(seq_sequence[i:i+60] + "\n")
                    seen_ids.add(current_seq_id)
                    combined_data[current_seq_id] = {
                        "name": current_seq_name,
                        "organism": current_seq_org,
                        "length": seq_length,
                        "families": [fam]
                    }
                else:
                    # Sekvensen er allerede lagt til, oppdater families-listen
                    if fam not in combined_data[current_seq_id]["families"]:
                        combined_data[current_seq_id]["families"].append(fam)
                count += 1
            # Parse header-linjen for ny sekvens
            header_line = line[1:]  # fjern '>'
            # UniProt FASTA header-format: sp|ACC|EntryName Protein Name OS=Organsime OX=TaxID ... 
            parts = header_line.split("|")
            if len(parts) >= 3:
                acc = parts[1]
                rest = parts[2]
            else:
                # Fallback hvis header ikke har forventet format
                header_parts = header_line.split()
                acc = header_parts[0]
                rest = " ".join(header_parts[1:]) if len(header_parts) > 1 else ""
            # Finn protein-navn og organisme i headeren
            os_idx = rest.find(" OS=")
            ox_idx = rest.find(" OX=")
            if os_idx != -1:
                name_section = rest[:os_idx]
            else:
                name_section = rest
            # EntryName (første ord) og proteinets fulle navn (resten av name_section)
            if " " in name_section:
                entry_name, protein_name = name_section.split(" ", 1)
            else:
                entry_name = name_section
                protein_name = ""
            organism = ""
            if os_idx != -1 and ox_idx != -1:
                organism = rest[os_idx+4 : ox_idx]  # tekst mellom 'OS=' og ' OX='
            # Start en ny sekvens
            current_seq_id = acc
            current_seq_name = protein_name
            current_seq_org = organism
            seq_buffer = []  # tøm bufferen for sekvensbokstaver
        else:
            # Dette er en linje med sekvensdata (aminosyresekvens)
            seq_buffer.append(line.strip())
    # Etter loop: avslutt siste sekvens for denne familien (hvis eksisterer)
    if current_seq_id is not None:
        seq_sequence = "".join(seq_buffer)
        seq_length = len(seq_sequence)
        # Skriv siste sekvens til familie-FASTA
        fam_file.write(f">{current_seq_id} {current_seq_name} [{current_seq_org}]\n")
        for i in range(0, len(seq_sequence), 60):
            fam_file.write(seq_sequence[i:i+60] + "\n")
        # Oppdater kombinert fil og metadata
        if current_seq_id not in seen_ids:
            combined_fasta_file.write(f">{current_seq_id} {current_seq_name} [{current_seq_org}]\n")
            for i in range(0, len(seq_sequence), 60):
                combined_fasta_file.write(seq_sequence[i:i+60] + "\n")
            seen_ids.add(current_seq_id)
            combined_data[current_seq_id] = {
                "name": current_seq_name,
                "organism": current_seq_org,
                "length": seq_length,
                "families": [fam]
            }
        else:
            if fam not in combined_data[current_seq_id]["families"]:
                combined_data[current_seq_id]["families"].append(fam)
        count += 1
    # Lukk familie-FASTA og loggfør antall
    fam_file.close()
    family_counts[fam] = count
    if count == 0:
        print(f"[ADVARSEL] Ingen sekvenser funnet for {fam} (sjekk familienavnet eller tilgjengelighet).")
    else:
        print(f"[INFO] Hentet {count} sekvenser for {fam}.")

# Lukk kombinert FASTA-fil
combined_fasta_file.close()

# Skriv ut metadata TSV-fil for alle unike sekvenser
metadata_filename = f"sequences_metadata_{date_str}.tsv"
with open(metadata_filename, "w") as meta_f:
    meta_f.write("SequenceID\tName\tOrganism\tLength\tFamilies\n")
    for seq_id, info in combined_data.items():
        name    = info["name"]
        org     = info["organism"]
        length  = info["length"]
        families_list = ", ".join(info["families"])
        meta_f.write(f"{seq_id}\t{name}\t{org}\t{length}\t{families_list}\n")

# Skriv ut JSON-fil med informasjon om kjøringen
run_metadata = {
    "run_id": run_id,
    "timestamp": timestamp_str,
    "families": families,
    "family_counts": family_counts,
    "total_sequences": sum(family_counts.values()),
    "unique_sequences": len(seen_ids),
    "data_sources": {
        "UniProt": "REST API (rest.uniprot.org) med CAZy-familiefilter",
        "CAZy": "Familieklassifisering fra cazy.org brukt i spørringer"
    }
}
json_filename = f"run_metadata_{date_str}.json"
with open(json_filename, "w") as json_f:
    json.dump(run_metadata, json_f, indent=4)

print(f"[SLUTT] Ferdig med kjøring {run_id}. Unike sekvenser: {len(seen_ids)} (totalt hentet: {sum(family_counts.values())}).")
print(f"[SLUTT] Se utdatafiler: {len(families)} familie-FASTA, 1 kombinert FASTA, 1 TSV, 1 JSON (datomerket {date_str}).")
