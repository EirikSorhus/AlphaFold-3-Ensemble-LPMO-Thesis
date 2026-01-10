import requests
import json

def debug_uniprot():
    url = "https://rest.uniprot.org/uniprotkb/search"
    query_id = "P12345"
    query = f"accession:{query_id}"
    
    # Original fields from the file
    fields = "accession,sequence,organism_name,protein_name,xref_interpro,ft_signal,ft_domain,ft_region,ft_binding,ft_site"
    
    params = {
        "query": query,
        "fields": fields,
        "format": "json",
        "size": 1
    }

    print(f"Testing UniProt API with ID: {query_id}")
    print(f"URL: {url}")
    print(f"Params: {json.dumps(params, indent=2)}")

    try:
        response = requests.get(url, params=params, timeout=30)
        
        print(f"\nStatus Code: {response.status_code}")
        print("Response Text:")
        print(response.text)
        
        if response.status_code == 200:
            print("\nSuccess! Parsed JSON:")
            try:
                print(json.dumps(response.json(), indent=2))
            except:
                print("Could not parse JSON.")
        else:
            print("\nRequest failed.")
            
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    debug_uniprot()
