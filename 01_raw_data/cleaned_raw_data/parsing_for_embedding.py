import json

with open('cleaned_uniprot_data.json', 'r') as file:
    data = json.load(file)

parsed_data = []

for entry in data.get("results", []):
    protein_id = entry.get("primaryAccession")
    protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
    sequence = entry.get("sequence", {}).get("value", "")
    
    go_annotations = []
    for cross_ref in entry.get("uniProtKBCrossReferences", []):
        if cross_ref.get("database") == "GO":
            for prop in cross_ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    go_annotations.append(prop.get("value"))

    parsed_data.append({
        "protein_id": protein_id,
        "protein_name": protein_name,
        "sequence": sequence,
        "go_annotations": go_annotations
    })

#print(json.dumps(parsed_data, indent=4))

with open('parsed_proteins.json', 'w') as outfile:
    json.dump(parsed_data, outfile, indent=4)

print(f"Data saved to parsed_proteins.json")