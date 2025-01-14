import json
import matplotlib.pyplot as plt
from collections import Counter


def analyze_go_terms(json_data):
    # Parse the JSON data (assuming it's already loaded as a Python dictionary)
    results = json_data['results']  # Assuming 'results' contains the list of annotations

    # List to store the GO terms
    go_terms = []

    # Iterate through the results to extract the GO terms
    for entry in results:
        # Ensure the entry has the dbReference type 'GO'
        if 'dbReference' in entry:
            for reference in entry['dbReference']:
                if reference.get('type') == 'GO':
                    # Extract the GO term from the 'id' field
                    go_terms.append(reference['id'])

    # Count the occurrences of each GO term
    go_term_counts = Counter(go_terms)

    # Sort the GO terms by frequency
    go_term_counts = dict(go_term_counts.most_common())

    # Plot the distribution
    plt.figure(figsize=(10, 6))
    plt.bar(go_term_counts.keys(), go_term_counts.values())
    plt.xlabel('GO Terms')
    plt.ylabel('Frequency')
    plt.title('Distribution of GO Terms')
    plt.xticks(rotation=90)  # Rotate x-axis labels for better visibility
    plt.tight_layout()  # Adjust layout to avoid clipping
    plt.show()

# Example usage:
# Assuming you have the JSON data loaded into `json_data`
# Example of loading JSON data (replace with actual JSON file or response):

json_file = ""
with open(json_file) as f:
    json_data = json.load(f)
    analyze_go_terms(json_data)
