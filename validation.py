from Bio.ExPASy import get_sprot_raw
from Bio.SwissProt import read
from config import DATA_DIR

def fetch_pfam_domains(uniprot_id):
    """
    Fetch Pfam domain annotations for a given UniProt ID using SwissProt.
    """
    try:
        # Fetch SwissProt record
        handle = get_sprot_raw(uniprot_id)
        record = read(handle)
        
        print(f"Pfam domains for {uniprot_id}:")
        for feature in record.features:
            if feature.type == "DOMAIN" and "Pfam" in feature.qualifiers.get("database", []):
                domain_name = feature.qualifiers.get("description", [""])[0].split(";")[0].strip()
                start = feature.location.start
                end = feature.location.end
                print(f"Domain: {domain_name}, Positions: {start}-{end}")
    except Exception as e:
        print(f"Error fetching Pfam domains: {e}")

def fetch_swissprot_annotations(uniprot_id):
    """
    Fetch SwissProt annotations for a given UniProt ID.
    """
    try:
        # Fetch SwissProt record
        handle = get_sprot_raw(uniprot_id)
        record = read(handle)
        
        print(f"SwissProt Annotations for {uniprot_id}:")
        
        # Extract function 
        if record.comments:
            for comment in record.comments:
                if comment.startswith("FUNCTION:"):
                    print(f"Function: {comment.replace('FUNCTION:', '').strip()}")
                    break
            else:
                print("Function: N/A")
        else:
            print("Function: N/A")
        
        # Extract subcellular location 
        if record.comments:
            for comment in record.comments:
                if comment.startswith("SUBCELLULAR LOCATION:"):
                    print(f"Subcellular Location: {comment.replace('SUBCELLULAR LOCATION:', '').strip()}")
                    break
            else:
                print("Subcellular Location: N/A")
        else:
            print("Subcellular Location: N/A")
    except Exception as e:
        print(f"Error fetching SwissProt annotations: {e}")


uniprot_id = "P38398"  # Human BRCA1
fetch_pfam_domains(uniprot_id)
fetch_swissprot_annotations(uniprot_id)