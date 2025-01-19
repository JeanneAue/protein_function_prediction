import json
import requests
from time import sleep

from obo_parsing import parse_obo_file


def get_hierarchy(go_terms, go_id):
    """
    Retrieve the complete hierarchy for a given GO term by its ID.

    Parameters:
        go_terms (dict): Dictionary of parsed GO terms.
        go_id (str): GO term ID to retrieve the hierarchy for.

    Returns:
        dict: A hierarchy dictionary for the given GO term.
    """
    if go_id not in go_terms:
        return None

    def traverse_hierarchy(term_id, visited):
        if term_id in visited:
            return {}
        visited.add(term_id)
        term = go_terms[term_id]
        hierarchy = {
            "id": term_id,
            "name": term.get("name", ""),
            "namespace": term.get("namespace", ""),
            "parents": [],
        }
        for parent_id in term.get("parents", []):
            hierarchy["parents"].append(traverse_hierarchy(parent_id, visited))
        return hierarchy

    return traverse_hierarchy(go_id, set())


def fetch_go_hierarchy_by_name(obo_file, go_name):
    """
    Fetch the hierarchy of a GO term by its name from an OBO file.

    Parameters:
        obo_file (str): Path to the OBO file.
        go_name (str): GO term name to retrieve the hierarchy for.

    Returns:
        dict: Hierarchy dictionary for the given GO term.
    """
    go_terms, name_to_id = parse_obo_file(obo_file)
    go_id = name_to_id.get(go_name)
    if not go_id:
        print(f"GO term with name '{go_name}' not found.")
        return None
    return get_hierarchy(go_terms, go_id)





def add_go_hierarchy_to_proteins(input_file, output_file, sleep_time=0.5):
    """
    Adds the GO term hierarchy to each protein in the JSON file.

    Parameters:
        input_file (str): Path to the input JSON file containing protein data.
        output_file (str): Path to the output JSON file with the GO hierarchy added.
        sleep_time (float): Time in seconds to wait between API calls (to avoid rate limiting).
    """
    obo_file = "./go-basic.obo"

    # Load the protein data
    with open(input_file, 'r') as f:
        proteins = json.load(f)

    for protein in proteins:
        # Check if the protein has GO annotations
        if "go_annotations" in protein and protein["go_annotations"]:
            protein["go_hierarchies"] = {}
            for go_annotation in protein["go_annotations"]:
                # Extract the GO term name (everything after the first colon)
                if ":" in go_annotation:
                    go_term_name = ":".join(go_annotation.split(":")[1:])  # Join all parts after the first colon
                    print(f"Fetching hierarchy for GO term: {go_term_name}")
                    hierarchy = fetch_go_hierarchy_by_name(obo_file, go_term_name)
                    if hierarchy:
                        protein.setdefault("go_hierarchies", {})[go_term_name] = hierarchy

    # Save the updated proteins to a new JSON file
    with open(output_file, 'w') as f:
        json.dump(proteins, f, indent=4)

    print(f"Updated data with GO hierarchies saved to {output_file}")

# Example usage
input_file = "./../raw_data/filtered_proteins.json"  # Replace with your input file path
output_file = "./../raw_data/proteins_with_complete_go_hierarchy.json"  # Replace with your desired output file path
add_go_hierarchy_to_proteins(input_file, output_file)

