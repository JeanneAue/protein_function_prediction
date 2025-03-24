import json
from obo_parsing import parse_obo_file
from prediction.debug import go_hierarchies


# not tested
def filter_proteins_by_parent_go_terms(input_file, output_file, obo_file, max_parents=6000):
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

    def get_parents(term_id, hierarchy):
        parents_of_term = set()
        if hierarchy:
            for element in hierarchy:
                if element:
                    if element["id"] == term_id:
                        for parent in element.get("parents", []):
                            if parent:
                                parent_id = parent.get("id")
                                parents_of_term.add(parent_id)
                    else:
                        parent_ids = get_parents(term_id, element["parents"])
                        parents_of_term.update(parent_ids)
        return parents_of_term

    def extract_go_parents_recursive(go_hierarchies, term_id, collected_parents=None):
        if collected_parents is None:
            collected_parents = set()

            # Add the current term to the collected parents
        collected_parents.add(term_id)

        # Get parent ids of term_id
        parent_ids = get_parents(term_id, go_hierarchies)

        for parent_id in parent_ids:
            if parent_id and parent_id not in collected_parents:
                extract_go_parents_recursive(go_hierarchies, parent_id, collected_parents)

        return collected_parents

    # Load the protein data
    with open(input_file, 'r') as f:
        proteins = json.load(f)

    filtered_proteins = []

    for protein in proteins:
        if "go_hierarchies" in protein and protein["go_hierarchies"]:
            go_hierarchy = protein["go_hierarchies"]
            all_parents = set()

            for go_term in go_hierarchy:
                go_id = go_term["id"]
                parent_terms = extract_go_parents_recursive(go_hierarchy, go_id)
                all_parents.update(parent_terms)

            # Add the protein to the filtered list if it has GO annotations
            if all_parents:
                filtered_proteins.append({
                    "protein_id": protein["protein_id"],
                    "protein_name": "no_name_now",
                    "sequence": "no_sequence",
                    "go_annotations": list(all_parents)
                })

    # Save the filtered proteins to a new JSON file
    with open(output_file, 'w') as f:
        json.dump(filtered_proteins, f, indent=4)

    print(f"Filtered proteins with GO annotations saved to {output_file}")


# Example usage
input_file = "../embeddings/protein_data_with_embeddings_and_hierarchy.json"  # Replace with your input file path
output_file = "../embeddings/compressed_protein_data.json"  # Replace with your desired output file path
obo_file = "./go-basic.obo"  # Path to the OBO file
max_parents = 5  # Specify the maximum number of parent GO terms to include

filter_proteins_by_parent_go_terms(input_file, output_file, obo_file)
