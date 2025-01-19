def count_lines_in_json(file_path):
    """
    Counts the number of lines in a JSON file.

    Parameters:
        file_path (str): Path to the JSON file.

    Returns:
        int: The total number of lines in the file.
    """
    try:
        with open(file_path, 'r') as file:
            line_count = sum(1 for _ in file)
        return line_count
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return 0
    except Exception as e:
        print(f"An error occurred: {e}")
        return 0


json_file_path = "./protein_data_with_embeddings_and_hierarchy.json"  # Replace with your JSON file path
lines = count_lines_in_json(json_file_path)

if lines > 0:
    print(f"The file '{json_file_path}' contains {lines} lines.")