from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re


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
        },
        "lengths": lengths
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
    fasta_file = "../01_raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.fasta"
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


def extract_data_from_stats(stats_dict):
    """Extract and prepare data from the statistics dictionary."""

    median = stats_dict['median_length']
    max_length = stats_dict['longest_sequence']['length']

    # lengths as ndarray
    lengths = np.array(stats_dict['lengths'])

    return {
        'lengths': lengths,
        'total': stats_dict['total_sequences'],
        'median': median,
        'min_length': stats_dict['shortest_sequence']['length'],
        'max_length': max_length,
        'shortest_id': stats_dict['shortest_sequence']['id'],
        'longest_id': stats_dict['longest_sequence']['id'],
        'shortest_desc': stats_dict['shortest_sequence']['description'],
        'longest_desc': stats_dict['longest_sequence']['description']
    }


def extract_protein_info(description):
    """Extract protein name and organism from the description."""
    match = re.search(r'(.*) OS=(.*) OX=', description)
    if match:
        protein = match.group(1)
        organism = match.group(2)
        return protein, organism
    return description, "Unknown"


def visualize_sequence_distribution(data):
    """Create a comprehensive visualization of sequence length distribution."""
    plt.figure(figsize=(14, 10))

    # 1. Logarithmic histogram of sequence lengths
    plt.subplot(2, 2, 1)

    # Use log scale for x-axis to better visualize the wide range
    bins = np.logspace(np.log10(max(1, data['min_length'])),
                       np.log10(data['max_length']),
                       50)

    sns.histplot(data['lengths'], bins=bins, kde=True, color='skyblue',
                 edgecolor='darkblue', alpha=0.7)

    plt.axvline(data['median'], color='red', linestyle='--',
                label=f"Median: {data['median']} aa")

    plt.xscale('log')
    plt.title("Distribution of Sequence Lengths (Log Scale)", fontsize=12)
    plt.xlabel("Sequence Length (amino acids)", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.legend()

    # 2. Sequence length categories
    plt.subplot(2, 2, 2)

    # Define length categories
    categories = [
        (0, 100, 'Very Short'),
        (100, 500, 'Short'),
        (500, 1000, 'Medium'),
        (1000, 5000, 'Long'),
        (5000, float('inf'), 'Very Long')
    ]

    # Count sequences in each category
    cat_counts = []
    cat_labels = []

    for min_len, max_len, label in categories:
        count = np.sum((data['lengths'] >= min_len) & (data['lengths'] < max_len))
        cat_counts.append(count)
        cat_labels.append(f"{label}\n({min_len}-{max_len if max_len != float('inf') else '∞'} aa)")

    # Plot
    plt.pie(cat_counts, labels=cat_labels, autopct='%1.1f%%',
            colors=sns.color_palette("coolwarm", len(categories)),
            startangle=90, explode=[0.05] * len(categories))
    plt.title("Sequence Length Categories", fontsize=12)

    # 3. Extreme sequences comparison
    plt.subplot(2, 2, 3)

    # Extract protein names for cleaner labels
    shortest_protein, shortest_organism = extract_protein_info(data['shortest_desc'])
    longest_protein, longest_organism = extract_protein_info(data['longest_desc'])

    # Format labels
    shortest_label = f"{shortest_protein.split('|')[-1]}\n({shortest_organism})"
    longest_label = f"{longest_protein.split('|')[-1]}\n({longest_organism})"

    # Create comparison bar chart with logarithmic scale
    plt.bar([shortest_label, longest_label],
            [data['min_length'], data['max_length']],
            color=['#8ecae6', '#e63946'], width=0.6)

    plt.yscale('log')
    plt.title("Shortest vs. Longest Sequence", fontsize=12)
    plt.ylabel("Sequence Length (log scale)", fontsize=10)

    # Add text annotations showing exact lengths
    plt.text(0, data['min_length'], f"{data['min_length']} aa",
             ha='center', va='bottom', fontsize=10)
    plt.text(1, data['max_length'] / 10, f"{data['max_length']} aa",
             ha='center', va='bottom', fontsize=10)

    # 4. Summary statistics
    plt.subplot(2, 2, 4)
    plt.axis('off')

    summary_text = (
        f"Dataset Summary:\n"
        f"------------------------\n"
        f"Total Sequences: {data['total']:,}\n"
        f"Median Length: {data['median']} amino acids\n"
        f"Shortest Sequence: {data['min_length']} amino acids\n"
        f"Longest Sequence: {data['max_length']} amino acids\n\n"
        f"Shortest: {shortest_protein.split('|')[-1]}\n"
        f"Longest: {longest_protein.split('|')[-1]} (Titin)\n\n"
        f"Length Range: {data['max_length'] / data['min_length']:,.1f}x difference"
    )

    plt.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=11,
             bbox=dict(facecolor='lavender', alpha=0.5, boxstyle='round,pad=1'))

    plt.tight_layout()
    return plt.gcf()


def compare_to_reference_databases(data):
    """Create a visualization comparing to reference protein databases."""
    plt.figure(figsize=(10, 6))

    # Reference database statistics (approximate values)
    databases = {
        'Your Dataset': data['median'],
        'SwissProt': 335,
        'TrEMBL': 317,
        'PDB': 250,
        'RefSeq': 375
    }

    # Create bar chart
    plt.bar(databases.keys(), databases.values(),
            color=sns.color_palette("muted", len(databases)))

    plt.title('Median Sequence Length Comparison\nwith Reference Databases', fontsize=14)
    plt.ylabel('Median Sequence Length (amino acids)', fontsize=12)
    plt.ylim(0, max(databases.values()) * 1.2)

    # Add values on top of bars
    for i, (db, value) in enumerate(databases.items()):
        plt.text(i, value + 20, str(value), ha='center', fontsize=10)

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    return plt.gcf()


# Main function to generate all visualizations
def generate_protein_sequence_visualizations(stats_dict):
    """Generate comprehensive visualizations for protein sequence data."""
    data = extract_data_from_stats(stats_dict)

    # Create and save all visualizations
    fig1 = visualize_sequence_distribution(data)
    fig2 = compare_to_reference_databases(data)

    return fig1, fig2


# Generate all visualizations
fig1, fig2 = generate_protein_sequence_visualizations(analysis)
plt.show()

# # Save the visualizations to files
fig1.savefig('sequence_distribution.png')