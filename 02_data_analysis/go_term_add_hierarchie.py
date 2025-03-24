import json
from obo_parsing import parse_obo_file


def get_hierarchy(go_terms, go_id, cached_hierarchies):
    """
    Retrieve the complete hierarchy for a given GO term by its ID.
    Caches the result for future use.

    Parameters:
        go_terms (dict): Dictionary of parsed GO terms.
        go_id (str): GO term ID to retrieve the hierarchy for.
        cached_hierarchies (dict): Cache for storing previously computed hierarchies.

    Returns:
        dict: A hierarchy dictionary for the given GO term.
    """
    if go_id in cached_hierarchies:
        return cached_hierarchies[go_id]

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

    hierarchy = traverse_hierarchy(go_id, set())
    cached_hierarchies[go_id] = hierarchy
    return hierarchy


def add_go_hierarchy_to_proteins(input_file, output_file, obo_file):
    """
    Adds the GO term hierarchy to each protein in the JSON file.

    Parameters:
        input_file (str): Path to the input JSON file containing protein data.
        output_file (str): Path to the output JSON file with the GO hierarchy added.
        obo_file (str): Path to the GO-basic OBO file.
    """
    # Load GO terms and create mapping
    go_terms, name_to_id = parse_obo_file(obo_file)

    # Cache for GO term hierarchies
    cached_hierarchies = {}

    # Load the protein data
    with open(input_file, 'r') as f:
        proteins = json.load(f)

    def build_hierarchy(go_id, visited):
        """Recursively build a consistent GO term hierarchy."""
        if go_id in visited:
            return None  # Avoid circular references
        visited.add(go_id)

        term = go_terms.get(go_id)
        if not term:
            return None

        return {
            "id": go_id,
            "name": term.get("name", ""),
            "namespace": term.get("namespace", ""),
            "parents": [
                build_hierarchy(parent_id, visited)
                for parent_id in term.get("parents", [])
                if parent_id in go_terms
            ]
        }

    for protein in proteins:
        if "go_annotations_with_ids" in protein and protein["go_annotations_with_ids"]:
            # Restructure hierarchies as a list of top-level terms
            protein["go_hierarchies"] = []
            for go_annotation in protein["go_annotations_with_ids"]:
                go_id = go_annotation["id"]

                if go_id:
                    print(f"Processing GO term: {go_id}")

                    # Check cache first to avoid recomputation
                    if go_id not in cached_hierarchies:
                        cached_hierarchies[go_id] = build_hierarchy(go_id, set())

                    hierarchy = cached_hierarchies[go_id]
                    if hierarchy:
                        protein["go_hierarchies"].append(hierarchy)

    with open(output_file, 'w') as f:
        json.dump(proteins, f, indent=4)

    print(f"Updated data with GO hierarchies saved to {output_file}")


# Example usage
input_file = "../embeddings/protein_embeddings_with_go_ids.json"
output_file = "../embeddings/protein_data_with_embeddings_and_hierarchy.json"
obo_file = "./go-basic.obo"

add_go_hierarchy_to_proteins(input_file, output_file, obo_file)
