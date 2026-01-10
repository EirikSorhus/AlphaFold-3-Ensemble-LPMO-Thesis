import requests
import time
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class InterProClient:
    """
    Client for fetching domain architecture information from the InterPro API.
    Specifically targets the /entry/interpro/protein/uniprot/{acc} endpoint
    to retrieve precise domain boundaries.
    """
    
    BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/"

    def __init__(self, max_retries=3, delay=1.0):
        self.max_retries = max_retries
        self.delay = delay

    def fetch_domains(self, uniprot_acc):
        """
        Fetches InterPro entries for a specific UniProt accession.
        Returns a list of dictionaries with domain details.
        """
        url = f"{self.BASE_URL}{uniprot_acc}"
        headers = {"Accept": "application/json"}
        
        attempt = 0
        while attempt < self.max_retries:
            try:
                response = requests.get(url, headers=headers, timeout=10)
                
                if response.status_code == 404:
                    # logger.warning(f"No InterPro data found for {uniprot_acc} (404).") # Common, reduce log spam
                    return []
                
                if response.status_code == 204:
                    # 204 No Content = Valid protein, but no InterPro entries
                    return []
                
                response.raise_for_status()
                data = response.json()
                
                return self._parse_response(data, uniprot_acc) 

            except json.JSONDecodeError:
                # Catch cases where 200 OK is returned but empty/bad body
                if response.status_code == 200 and not response.text.strip():
                     return []
                logger.error(f"Invalid JSON/Empty response for {uniprot_acc}. Status: {response.status_code}")
                return []


            except requests.exceptions.RequestException as e:
                attempt += 1
                logger.warning(f"Attempt {attempt}/{self.max_retries} failed for {uniprot_acc}: {e}")
                time.sleep(self.delay * attempt)
        
        logger.error(f"Failed to fetch InterPro data for {uniprot_acc} after {self.max_retries} retries.")
        return []

    def _parse_response(self, data, target_acc):
        """
        Parses the raw JSON response from InterPro into a simplified list of features.
        """
        parsed_entries = []
        
        # Results is usually a list of InterPro entries (families, domains, etc.) that match the protein
        results = data.get("results", [])
        
        for res in results:
            metadata = res.get("metadata", {})
            ipr_acc = metadata.get("accession")
            ipr_name = metadata.get("name")
            ipr_type = metadata.get("type") # e.g. Domain, Family, Homologous_superfamily
            
            # The 'proteins' list contains matches. Since we queried BY protein, 
            # we expect our target_acc to be in there.
            proteins = res.get("proteins", [])
            
            for prot in proteins:
                # API is case-insensitive but returns lower case usually
                if prot.get("accession").lower() == target_acc.lower():
                    
                    locations = prot.get("entry_protein_locations", [])
                    for loc in locations:
                        fragments = loc.get("fragments", [])
                        for frag in fragments:
                            parsed_entries.append({
                                "id": ipr_acc,
                                "name": ipr_name,
                                "type": ipr_type,
                                "start": frag.get("start"),
                                "end": frag.get("end"),
                                "status": frag.get("dc-status") # e.g. CONTINUOUS
                            })
                            
        return parsed_entries

if __name__ == "__main__":
    # Internal test
    client = InterProClient()
    print("Testing InterProClient with W4FSN3...")
    domains = client.fetch_domains("W4FSN3")
    for d in domains:
        print(d)
