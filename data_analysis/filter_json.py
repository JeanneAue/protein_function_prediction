import json

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

# Example usage
input_file = "./../raw_data/parsed_proteins.json"  # Replace with your input file path
output_file = "./../raw_data/filtered_proteins.json"  # Replace with your desired output file path
filter_proteins_with_go_terms(input_file, output_file)
