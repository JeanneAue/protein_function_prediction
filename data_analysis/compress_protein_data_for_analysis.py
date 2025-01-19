import json
from obo_parsing import parse_obo_file

# not tested
def filter_proteins_by_parent_go_terms(input_file, output_file, obo_file, max_parents):
    """
    Creates a new JSON file including only a specified number of parent GO terms
    for each protein. The output includes only GO term IDs in the "go_annotations" field.

    Parameters:
        input_file (str): Path to the input JSON file containing protein data.
        output_file (str): Path to the output JSON file with filtered GO annotations.
        obo_file (str): Path to the GO-basic OBO file.
        max_parents (int): Maximum number of parent GO terms to include for each annotation.
    """
    # Load GO terms and create mapping
    go_terms, _ = parse_obo_file(obo_file)

    def get_parent_terms(go_id, max_depth):
        """
        Retrieve parent terms up to a specified depth for a given GO term ID.

        Parameters:
            go_id (str): The GO term ID.
            max_depth (int): Maximum number of parent terms to include.

        Returns:
            list: List of parent GO term IDs up to the specified depth.
        """
        if go_id not in go_terms:
            return []

        parent_terms = []
        to_visit = [(go_id, 0)]  # (current GO ID, current depth)
        visited = set()

        while to_visit:
            current_id, depth = to_visit.pop(0)
            if current_id in visited or depth >= max_depth:
                continue
            visited.add(current_id)
            parent_terms.append(current_id)

            # Add parents of the current term to the visit list
            parents = go_terms[current_id].get("parents", [])
            to_visit.extend((parent_id, depth + 1) for parent_id in parents)

        return parent_terms

    # Load the protein data
    with open(input_file, 'r') as f:
        proteins = json.load(f)

    filtered_proteins = []

    for protein in proteins:
        if "go_annotations" in protein and protein["go_annotations"]:
            new_go_annotations = set()
            for go_annotation in protein["go_annotations"]:
                if ":" in go_annotation:
                    go_term_name = ":".join(go_annotation.split(":")[1:])
                    go_id = go_terms.get(go_term_name, {}).get("id")
                    if go_id:
                        # Get parent terms and limit to max_parents
                        parent_terms = get_parent_terms(go_id, max_parents)
                        new_go_annotations.update(parent_terms)

            # Add the protein to the filtered list if it has GO annotations
            if new_go_annotations:
                filtered_proteins.append({
                    "protein_id": protein["protein_id"],
                    "go_annotations": list(new_go_annotations)
                })

    # Save the filtered proteins to a new JSON file
    with open(output_file, 'w') as f:
        json.dump(filtered_proteins, f, indent=4)

    print(f"Filtered proteins with GO annotations saved to {output_file}")


# Example usage
input_file = "../embeddings/protein_data_with_embeddings_and_hierarchie.json"  # Replace with your input file path
output_file = "../embeddings/compressed_protein_data.json"  # Replace with your desired output file path
obo_file = "./go-basic.obo"  # Path to the OBO file
max_parents = 5  # Specify the maximum number of parent GO terms to include

filter_proteins_by_parent_go_terms(input_file, output_file, obo_file, max_parents)
