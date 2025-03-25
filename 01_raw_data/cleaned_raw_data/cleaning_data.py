import json

def clean_uniprot_json(input_json):
    """
    Cleans up UniProt JSON data.
    """

    data = json.loads(json.dumps(input_json))


    results = data.get("results", [])

    for entry in results:
        for key_to_remove in ["entryAudit", "organism", "references"]:
            if key_to_remove in entry:
                del entry[key_to_remove]

        # 2) Filter uniProtKBCrossReferences
        if "uniProtKBCrossReferences" in entry:
            filtered_cross_refs = []
            for cross_ref in entry["uniProtKBCrossReferences"]:
                
                props = cross_ref.get("properties", [])
                keep_this = any(
                    p.get("key") == "GoTerm" and p.get("value", "").startswith("F:")
                    for p in props
                )
                if keep_this:
                    filtered_cross_refs.append(cross_ref)

            # Overwrite the original list
            entry["uniProtKBCrossReferences"] = filtered_cross_refs

    return data


if __name__ == "__main__":
    file_name = '01_raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.json'
    with open(file_name, 'r', encoding="utf-8") as f1:
            input_data = json.load(f1)

    cleaned_data = clean_uniprot_json(input_data)

    output_file = "01_raw_data/cleaned_raw_data/cleaned_uniprot_data.json"
    with open(output_file, "w", encoding="utf-8") as outfile:
        json.dump(cleaned_data, outfile, indent=2)

    print(f"Cleaned JSON data saved to {output_file}")
