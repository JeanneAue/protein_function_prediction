import json

from obo_parsing import name_to_id

def filter_proteins_with_go_terms(input_file, output_file):
    """
    Filters proteins with at least one GO term and writes them to a new JSON file.

    Parameters:
        input_file (str): Path to the input JSON file containing protein data.
        output_file (str): Path to the output JSON file to save the filtered data.
    """
    # Load the JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)

    # Filter proteins that have at least one GO term
    filtered_proteins = [protein for protein in data if protein.get("go_annotations")]

    # Write the filtered proteins to a new JSON file
    with open(output_file, 'w') as f:
        json.dump(filtered_proteins, f, indent=4)

    print(f"Filtered data saved to {output_file} with {len(filtered_proteins)} proteins.")

def filter_proteins_with_go_ids(input_file, output_file, name_to_id):
    """
    Filter proteins to only retain 'go_annotations_with_ids' field based on GO term names.

    Parameters:
        input_file (str): Path to the input JSON file containing protein data.
        output_file (str): Path to save the filtered JSON data.
        name_to_id (dict): Dictionary mapping GO term names to their corresponding IDs.
    """
    with open(input_file, 'r') as f:
        data = json.load(f)

    filtered_data = []

    for protein in data:
        # Prepare a new protein entry with only protein_id and go_annotations_with_ids
        filtered_protein = {
            "protein_id": protein["protein_id"],
            "protein_name": protein["protein_name"],
            "sequence": protein["sequence"],
            "go_annotations_with_ids": []
        }

        for annotation in protein.get("go_annotations", []):
            # Extract the GO term name (everything after the first colon)
            if ":" in annotation:
                go_term_name = ":".join(annotation.split(":")[1:])

                # Find the corresponding GO ID
                go_id = name_to_id.get(go_term_name)

                # Append the annotation and its ID as a dictionary if ID exists
                if go_id:
                    filtered_protein["go_annotations_with_ids"].append({
                        "name": go_term_name,
                        "id": go_id
                    })

        # Add protein only if it has at least one valid GO annotation with an ID
        if filtered_protein["go_annotations_with_ids"]:
            filtered_data.append(filtered_protein)

    with open(output_file, 'w') as f:
        json.dump(filtered_data, f, indent=4)

    print(f"Filtered data saved to {output_file} with {len(filtered_data)} proteins.")


# Example usage
input_file = "../01_raw_data/cleaned_raw_data/parsed_proteins.json"  # Replace with your input file path
output_file = "../01_raw_data/cleaned_raw_data/filtered_proteins.json"
output_file_2 = "../01_raw_data/cleaned_raw_data/filtered_proteins_with_GO_ids.json"
filter_proteins_with_go_ids(output_file, output_file_2, name_to_id)
