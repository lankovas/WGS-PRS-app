import os
import requests
import gzip

# ==============================================================================
# PGS CATALOG DOWNLOAD UTILITY
# ==============================================================================

def fetch_pgs_panel(pgs_id, download_dir="../Pipeline/inputs"):
    """
    Downloads a GRCh38 harmonized Polygenic Risk Score (PRS) panel from the PGS Catalog.
    Also retrieves the associated metadata (trait, variant count, and ancestry distribution)
    directly from the REST API to provide context in the UI.
    """
    # Ensure the target directory exists before downloading
    os.makedirs(download_dir, exist_ok=True)
    
    file_name = f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    file_path = os.path.join(download_dir, file_name)
    
    # Endpoints for downloading the actual file and fetching the metadata
    ftp_url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{file_name}"
    api_url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    
    try:
        # ----------------------------------------------------------------------
        # 1. Fetch Metadata from the REST API
        # ----------------------------------------------------------------------
        api_response = requests.get(api_url)
        
        # Initialize default information in case the API call fails or is incomplete
        info = {
            "trait": "Unknown", 
            "variants": "Unknown", 
            "gwas_ancestry": "Not reported",
            "training_ancestry": "Not reported",
            "evaluation_ancestry": "Not reported",
            "pgs_id": pgs_id
        }
        
        if api_response.status_code == 200:
            api_data = api_response.json()
            
            # The API provides the exact ancestry percentages shown on the website
            ancestry_dist = api_data.get("ancestry_distribution", {})
            
            # Map standard 3-letter ancestry codes to readable labels
            anc_map = {
                "EUR": "European", "EAS": "East Asian", "AFR": "African", 
                "SAS": "South Asian", "AMR": "Hispanic/Latin American", 
                "GME": "Greater Middle Eastern", "OTH": "Other", 
                "NR": "Not Reported", "MAC": "Multi-Ancestry"
            }
            
            def parse_distribution(category_key):
                """Helper function to format the ancestry percentages for a specific category (e.g., GWAS)."""
                dist_dict = ancestry_dist.get(category_key, {}).get("dist", {})
                if not dist_dict:
                    return "Not reported for this score"
                
                # Sort the ancestries by percentage in descending order
                sorted_dist = sorted(dist_dict.items(), key=lambda item: item[1], reverse=True)
                
                formatted = []
                for k, v in sorted_dist:
                    name = anc_map.get(k, k) 
                    # Format string to 1 decimal place, e.g., "European (100.0%)"
                    formatted.append(f"{name} ({float(v):.1f}%)")
                    
                return ", ".join(formatted)

            # Apply the parser to the three main stages of panel development
            info["gwas_ancestry"] = parse_distribution("gwas")
            info["training_ancestry"] = parse_distribution("training")
            info["evaluation_ancestry"] = parse_distribution("eval")
            
        # ----------------------------------------------------------------------
        # 2. Download the Harmonized File
        # ----------------------------------------------------------------------
        ftp_response = requests.get(ftp_url, stream=True)
        
        if ftp_response.status_code == 404:
            return False, f"ID '{pgs_id}' not found, or no GRCh38 harmonized file is available."
        elif ftp_response.status_code != 200:
            return False, f"Failed to connect to PGS Catalog (Status: {ftp_response.status_code})"
            
        # Stream the file to disk in chunks to avoid memory issues with large panels
        with open(file_path, "wb") as f:
            for chunk in ftp_response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        # ----------------------------------------------------------------------
        # 3. Parse Local File Header
        # ----------------------------------------------------------------------
        # Read the top of the newly downloaded file to extract trait and variant count
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                if line.startswith('#trait_reported='):
                    info['trait'] = line.strip().split('=')[1]
                elif line.startswith('#variants_number='):
                    info['variants'] = line.strip().split('=')[1]
                    
        return True, {"file_path": file_path, "info": info}
        
    except Exception as e:
        return False, f"An error occurred: {str(e)}"