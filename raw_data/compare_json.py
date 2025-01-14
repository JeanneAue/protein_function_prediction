import json

def compare_json_files(file1, file2, char_limit=100000000):
    try:
        # Load JSON files
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            json1 = json.load(f1)
            json2 = json.load(f2)
        
        # Convert JSON data to strings
        json_str1 = json.dumps(json1, sort_keys=True, indent=4)
        json_str2 = json.dumps(json2, sort_keys=True, indent=4)

        # Truncate strings to the specified character limit
        json_str1 = json_str1[:char_limit]
        json_str2 = json_str2[:char_limit]

        # Compare characters
        non_matching = []
        for i, (char1, char2) in enumerate(zip(json_str1, json_str2)):
            if char1 != char2:
                non_matching.append((i, char1, char2))

        # Calculate additional non-matching characters if lengths differ
        if len(json_str1) != len(json_str2):
            longer, shorter = (json_str1, json_str2) if len(json_str1) > len(json_str2) else (json_str2, json_str1)
            extra_chars = longer[len(shorter):]
            non_matching.extend([(i, char, None) for i, char in enumerate(extra_chars, start=len(shorter))])

        # Total non-matching characters
        total_non_matching = len(non_matching)

        # Show first 300 non-matching characters
        first_300_mismatches = non_matching[:300]

        # Output results
        print(f"Total non-matching characters (limited to {char_limit} characters): {total_non_matching}")
        print("First 300 non-matching characters (index, file1_char, file2_char):")
        for mismatch in first_300_mismatches:
            print(mismatch)

    except Exception as e:
        print(f"Error: {e}")

# Usage example
# Replace 'file1.json' and 'file2.json' with the paths to your JSON files
compare_json_files('raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.json', 'raw_data/uniprotkb_proteome_UP000005640_AND_revi_2025_01_14.json')
