from Bio import SeqIO

fasta_file = "../raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.fasta"

# Parse the FASTA file
# each record parsed by SeqIO has the following attributes
#   - annotation (always empty)
#   - description (basically entire entry execpt for the sequence
#   - features
content = SeqIO.parse(fasta_file, "fasta")
count = 0
for record in content:
    count += 1
    if record.features:
        print(record.features)
    # print(f"ID: {record.id}")
    # print(f"Description: {record.description}")
    # print(f"Sequence: {record.seq[:50]}...")
    # print(f"Sequence length: {len(record.seq)}")
print(f"Entire length: {count}")