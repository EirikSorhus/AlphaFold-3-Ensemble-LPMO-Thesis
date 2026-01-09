import re
import requests
from typing import Dict, Optional, Tuple


def parse_ncbi_header(header_line: str) -> Dict[str, str]:
	line = header_line.lstrip(">")
	accession = line.split()[0]
	organism = ""
	org_match = re.search(r"\[(.*?)\]\s*$", line)
	if org_match:
		organism = org_match.group(1)
	protein_name = line.replace(accession, "").replace(f"[{organism}]", "").strip()
	return {
		"UniProt_ID": accession,
		"Protein_Name": protein_name,
		"Gene_Name": "NA",
		"Organism": organism,
		"EC_Number": "NA",
		"Reviewed": "unreviewed",
		"PDB_IDs": "",
		"Pfam_IDs": "",
		"InterPro_IDs": "",
		"Match_Status": "matched_ncbi",
		"Source_DB": "NCBI_Fallback",
	}


def fetch_from_ncbi_with_metadata(accession: str, contact_email: Optional[str] = None) -> Tuple[Optional[str], Optional[Dict[str, str]]]:
	base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
	params = {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}
	if contact_email:
		params["email"] = contact_email

	try:
		r = requests.get(base_url, params=params, timeout=20)
		if r.ok and ">" in r.text:
			lines = r.text.splitlines()
			header = lines[0]
			sequence = "".join(lines[1:])
			metadata = parse_ncbi_header(header)
			return sequence, metadata
	except Exception as exc:
		print(f"      [!] NCBI fetch failed for {accession}: {exc}")
	return None, None
