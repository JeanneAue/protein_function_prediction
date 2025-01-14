from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import torch



def analyze_and_visualize_fasta(file_path):
    """
    Analyze a FASTA file to determine the median, longest, and shortest sequences,
    and visualize the results.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with statistics about the sequences.
    """
    # Parse the sequences from the FASTA file
    sequences = list(SeqIO.parse(file_path, "fasta"))

    if not sequences:
        print("No sequences found in the FASTA file.")
        return None

    # Calculate lengths of sequences
    lengths = [len(record.seq) for record in sequences]

    # Identify the shortest and longest sequences
    shortest_idx = lengths.index(min(lengths))
    longest_idx = lengths.index(max(lengths))

    # Calculate the median length
    median_length = int(np.median(lengths))

    # Prepare results
    result = {
        "total_sequences": len(sequences),
        "median_length": median_length,
        "shortest_sequence": {
            "id": sequences[shortest_idx].id,
            "description": sequences[shortest_idx].description,
            "length": lengths[shortest_idx],
            "sequence": str(sequences[shortest_idx].seq),
        },
        "longest_sequence": {
            "id": sequences[longest_idx].id,
            "description": sequences[longest_idx].description,
            "length": lengths[longest_idx],
            "sequence": str(sequences[longest_idx].seq),
        }
    }

    # Visualize the results
    visualize_sequence_lengths(lengths, result)

    return result


def visualize_sequence_lengths(lengths, stats):
    """
    Visualize sequence length statistics.

    Args:
        lengths (list): List of sequence lengths.
        stats (dict): Statistics dictionary with shortest and longest sequences.
    """
    # Plot histogram of sequence lengths
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.hist(lengths, bins=30, color='skyblue', edgecolor='black')
    plt.axvline(stats["median_length"], color='red', linestyle='--', label="Median")
    plt.title("Distribution of Sequence Lengths")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.legend()

    # Plot bar chart for shortest and longest sequences
    plt.subplot(1, 2, 2)
    ids = [stats["shortest_sequence"]["id"], stats["longest_sequence"]["id"]]
    lengths = [stats["shortest_sequence"]["length"], stats["longest_sequence"]["length"]]
    colors = ['green', 'orange']
    plt.bar(ids, lengths, color=colors)
    plt.title("Shortest and Longest Sequences")
    plt.ylabel("Sequence Length")
    plt.xlabel("Sequence ID")

    # Show plots
    plt.tight_layout()
    plt.show()


# Example usage
if __name__ == "__main__":
    fasta_file = "../raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.fasta"
    analysis = analyze_and_visualize_fasta(fasta_file)
    if analysis:
        print(analysis)


# # Parse the FASTA file
# content = SeqIO.parse(fasta_file, "fasta")
# sequences = []
# for record in content:
#     print(f"ID: {record.id}")
#     print(f"Description: {record.description}")
#     print(f"Sequence: {record.seq[:50]}...")
#     sequences.append(record.seq)
#     print(f"Sequence length: {len(record.seq)}")
#
# print(f"Total length: {len(sequences)}")
# print(f"Median length of sequence: {np.median(sequences)}")