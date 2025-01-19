

def parse_obo_file(obo_file):
    """
    Parse a GO OBO file to extract term relationships and build a name-to-ID mapping.

    Parameters:
        obo_file (str): Path to the OBO file.

    Returns:
        dict: A dictionary of GO terms keyed by their IDs.
        dict: A reverse mapping from GO term names to their IDs.
    """
    go_terms = {}
    name_to_id = {}
    current_term = None

    with open(obo_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current_term = {}
            elif not line and current_term:
                if "id" in current_term:
                    go_terms[current_term["id"]] = current_term
                    if "name" in current_term:
                        name_to_id[current_term["name"]] = current_term["id"]
                current_term = None
            elif current_term is not None:
                if line.startswith("id: "):
                    current_term["id"] = line[4:]
                elif line.startswith("name: "):
                    current_term["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current_term["namespace"] = line[11:]
                elif line.startswith("def: "):
                    current_term["def"] = line[5:]
                elif line.startswith("is_a: "):
                    parent_id = line.split()[1]
                    current_term.setdefault("parents", []).append(parent_id)

    return go_terms, name_to_id

obo_file = "./go-basic.obo"
go_terms, name_to_id = parse_obo_file(obo_file)

#print(name_to_id)
print(name_to_id['fructose:sodium symporter activity'])
