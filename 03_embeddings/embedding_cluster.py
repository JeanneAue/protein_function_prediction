from transformers import BertTokenizer, BertModel
from torch.cuda.amp import autocast
import torch
import numpy as np
import json

tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
model = BertModel.from_pretrained("Rostlab/prot_bert")

if torch.cuda.is_available():
    model.to("cuda")

def generate_embedding(sequences, batch_size=8):
    all_embeddings = []
    global counter
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i : i + batch_size]
        spaced_sequences = [" ".join(seq) for seq in batch]

        inputs = tokenizer(spaced_sequences, return_tensors="pt", padding=True, truncation=True, max_length=1024)

        counter += 1
        print(f"Repetition: {counter}")
        print(f"Percentage: {counter / (len(sequences) // batch_size) * 100:.2f}%")
        print("--------------")

        if torch.cuda.is_available():
            inputs = {key: val.to("cuda") for key, val in inputs.items()}

        with torch.no_grad():
            with autocast():
                outputs = model(**inputs)
                batch_embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()

        all_embeddings.append(batch_embeddings)

    return np.vstack(all_embeddings)

with open("parsed_proteins.json", "r") as infile:
    parsed_data = json.load(infile)

sequences = [entry["sequence"] for entry in parsed_data]
protein_ids = [entry["protein_id"] for entry in parsed_data]
go_annotations_with_ids = [entry["go_annotations_with_ids"] for entry in parsed_data]

batch_size = 8
counter = 0
print("Generating embeddings for protein sequences...")
embeddings = generate_embedding(sequences, batch_size=batch_size)

protein_embeddings = []
for i, embedding in enumerate(embeddings):
    go_terms = [{"name": go["name"], "id": go["id"]} for go in go_annotations_with_ids[i]]
    protein_embeddings.append({
        "protein_id": protein_ids[i],
        "embedding": embedding.tolist(),
        "go_annotations_with_ids": go_terms
    })

with open("protein_embeddings_with_go_ids.json", "w") as outfile:
    json.dump(protein_embeddings, outfile, indent=4)

print("Protein embeddings saved to 'protein_embeddings_with_go_ids.json'")
