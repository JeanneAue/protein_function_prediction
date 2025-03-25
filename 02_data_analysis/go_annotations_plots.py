import json
from collections import defaultdict, Counter
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def GO_term_distribution_top_n(json_data, top_n=20):
    """
    Analyze and plot the distribution of GO terms from JSON data.

    Parameters:
        json_data (dict): Parsed JSON data containing GO annotations.
        top_n (int): Number of top GO terms to display in the plot. Default is 20.
    """
    # Parse the JSON data
    results = json_data['results']  # Assuming 'results' contains the list of annotations

    # List to store the GO terms
    go_terms = []

    # Iterate through the results to extract the GO terms
    for entry in results:
        if 'uniProtKBCrossReferences' in entry:
            for reference in entry['uniProtKBCrossReferences']:
                if 'properties' in reference:
                    for prop in reference['properties']:
                        if prop['key'] == "GoTerm":
                            go_terms.append(reference['id'][:15])

    # Count the occurrences of each GO term
    go_term_counts = Counter(go_terms)

    # Assuming that if more than top 100 are wanted, default to show all
    if top_n > 100:
        top_n = len(go_terms)

    # Sort and select the top N GO terms
    go_term_counts = dict(go_term_counts.most_common(top_n))

    # Plot the distribution as a horizontal bar chart
    plt.figure(figsize=(12, 8))  # Larger figure size
    bars = plt.barh(list(go_term_counts.keys()), list(go_term_counts.values()), color='skyblue')
    if top_n > 50:
        plt.gca().set_yticks([])

    plt.xlabel('Frequency')
    plt.ylabel('GO Terms')
    plt.title(f'Distribution of Top {top_n} GO Terms')

    if top_n <= 50:
        # Add labels next to the end of each bar
        for bar in bars:
            plt.text(
                bar.get_width(),  # x-coordinate: end of the bar
                bar.get_y() + bar.get_height() / 2,  # y-coordinate: center of the bar
                f'{int(bar.get_width())}',  # Text: width of the bar
                va='center',  # Vertically center-align the text
                ha='left'  # Horizontally align text to the left of the bar end
            )

    plt.tight_layout()  # Adjust layout to avoid clipping
    plt.gca().invert_yaxis()  # Invert y-axis for better alignment
    plt.show()

def GO_term_cooccurrence_heatmap(json_data, top_n=20):
    """
    Analyze and visualize GO term co-occurrence from JSON data.

    Parameters:
        json_data (dict): Parsed JSON data containing GO annotations.
        top_n (int): Number of top GO terms to analyze and display. Default is 20.
    """
    # Parse the JSON data
    results = json_data['results']  # Assuming 'results' contains the list of annotations

    # Dictionary to count co-occurrences of GO terms
    go_term_pairs = defaultdict(int)
    go_terms_per_entry = []

    # Iterate through the results to extract the GO terms
    for entry in results:
        go_terms = set()
        if 'uniProtKBCrossReferences' in entry:
            for reference in entry['uniProtKBCrossReferences']:
                if 'properties' in reference:
                    for prop in reference['properties']:
                        if prop['key'] == "GoTerm":
                            go_terms.add(reference['id'][:15])  # Use truncated ID
        if go_terms:
            go_terms_per_entry.append(go_terms)

    # Count co-occurrences for all unique pairs of GO terms
    for go_terms in go_terms_per_entry:
        for pair in itertools.combinations(go_terms, 2):
            go_term_pairs[tuple(sorted(pair))] += 1

    # Convert co-occurrence counts to a DataFrame for visualization
    cooccurrence_df = pd.DataFrame.from_dict(go_term_pairs, orient="index", columns=["count"])
    cooccurrence_df.reset_index(inplace=True)
    cooccurrence_df[["GO1", "GO2"]] = pd.DataFrame(cooccurrence_df["index"].tolist(), index=cooccurrence_df.index)
    cooccurrence_df = cooccurrence_df.drop(columns=["index"])

    # Filter to top N most frequent GO terms
    all_counts = defaultdict(int)
    for count, go1, go2 in cooccurrence_df.values:
        all_counts[go1] += count
        all_counts[go2] += count
    top_go_terms = [go for go, _ in Counter(all_counts).most_common(top_n)]
    filtered_df = cooccurrence_df[
        cooccurrence_df["GO1"].isin(top_go_terms) & cooccurrence_df["GO2"].isin(top_go_terms)
    ]

    # Pivot into a matrix for heatmap
    heatmap_df = filtered_df.pivot(index="GO1", columns="GO2", values="count").fillna(0)
    heatmap_df += heatmap_df.T  # Symmetrize the matrix

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_df, annot=False, fmt=".0f", cmap="Blues", square=True, cbar_kws={"label": "Co-occurrence Count"})
    plt.title(f"GO Term Co-occurrence Heatmap (Top {top_n} Terms)", fontsize=14)
    plt.xlabel("GO Terms")
    plt.ylabel("GO Terms")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

def GO_term_cooccurrence_heatmap_with_diagonal(json_data, top_n=20):
    """
    Analyze and visualize GO term co-occurrence from JSON data, including self-co-occurrences.

    Parameters:
        json_data (dict): Parsed JSON data containing GO annotations.
        top_n (int): Number of top GO terms to analyze and display. Default is 20.
    """
    # Parse the JSON data
    results = json_data['results']  # Assuming 'results' contains the list of annotations

    # Dictionary to count co-occurrences of GO terms
    go_term_pairs = defaultdict(int)
    go_terms_per_entry = []

    # Iterate through the results to extract the GO terms
    for entry in results:
        go_terms = set()
        if 'uniProtKBCrossReferences' in entry:
            for reference in entry['uniProtKBCrossReferences']:
                if 'properties' in reference:
                    for prop in reference['properties']:
                        if prop['key'] == "GoTerm":
                            go_terms.add(reference['id'][:15])  # Use truncated ID
        if go_terms:
            go_terms_per_entry.append(go_terms)

    # Count co-occurrences for all unique pairs of GO terms (including self-cooccurrence)
    for go_terms in go_terms_per_entry:
        for pair in itertools.combinations_with_replacement(go_terms, 2):  # Allow same-term pairs
            go_term_pairs[tuple(sorted(pair))] += 1

    # Convert co-occurrence counts to a DataFrame for visualization
    cooccurrence_df = pd.DataFrame.from_dict(go_term_pairs, orient="index", columns=["count"])
    cooccurrence_df.reset_index(inplace=True)
    cooccurrence_df[["GO1", "GO2"]] = pd.DataFrame(cooccurrence_df["index"].tolist(), index=cooccurrence_df.index)
    cooccurrence_df = cooccurrence_df.drop(columns=["index"])

    # Filter to top N most frequent GO terms
    all_counts = defaultdict(int)
    for count, go1, go2 in cooccurrence_df.values:
        all_counts[go1] += count
        all_counts[go2] += count
    top_go_terms = [go for go, _ in Counter(all_counts).most_common(top_n)]
    filtered_df = cooccurrence_df[
        cooccurrence_df["GO1"].isin(top_go_terms) & cooccurrence_df["GO2"].isin(top_go_terms)
    ]

    # Pivot into a matrix for heatmap
    heatmap_df = filtered_df.pivot(index="GO1", columns="GO2", values="count").fillna(0)
    heatmap_df += heatmap_df.T  # Symmetrize the matrix
    for term in heatmap_df.columns:  # Add diagonal counts (self-cooccurrence)
        heatmap_df.loc[term, term] /= 2  # Undo double-counting of diagonal elements

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_df, annot=False, fmt=".0f", cmap="Blues", square=True, cbar_kws={"label": "Co-occurrence Count"})
    plt.title(f"GO Term Co-occurrence Heatmap (Top {top_n} Terms)", fontsize=14)
    plt.xlabel("GO Terms")
    plt.ylabel("GO Terms")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def improved_GO_term_cooccurrence_heatmap(json_data, top_n=20, include_diagonal=False,
                                          log_scale=True, include_annotations=True,
                                          cluster_terms=True, figsize=(14, 12)):
    """
    Analyze and visualize GO term co-occurrence from JSON data with improved visualization.

    Parameters:
        json_data (dict): Parsed JSON data containing GO annotations.
        top_n (int): Number of top GO terms to analyze and display. Default is 20.
        include_diagonal (bool): Whether to include self-co-occurrences. Default is False.
        log_scale (bool): Whether to use log scale for color mapping. Default is True.
        include_annotations (bool): Whether to show count values in cells. Default is True.
        cluster_terms (bool): Whether to cluster similar terms together. Default is True.
        figsize (tuple): Figure size as (width, height). Default is (14, 12).
    """
    import json
    from collections import defaultdict, Counter
    import itertools
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from scipy.cluster import hierarchy

    # Parse the JSON data
    results = json_data['results']  # Assuming 'results' contains the list of annotations

    # Dictionary to store GO term descriptions
    go_term_descriptions = {}

    # Dictionary to count co-occurrences of GO terms
    go_term_pairs = defaultdict(int)
    go_terms_per_entry = []

    # Iterate through the results to extract the GO terms
    for entry in results:
        go_terms = set()
        if 'uniProtKBCrossReferences' in entry:
            for reference in entry['uniProtKBCrossReferences']:
                if reference.get('database') == 'GO':  # Make sure it's a GO reference
                    go_id = reference['id'][:15]  # Use truncated ID
                    go_terms.add(go_id)

                    # Extract GO term description if available
                    if 'properties' in reference:
                        for prop in reference['properties']:
                            if prop['key'] == "GoTerm":
                                go_term_descriptions[go_id] = prop['value']

        if go_terms:
            go_terms_per_entry.append(go_terms)

    # Count occurrences of individual GO terms
    go_term_counts = Counter()
    for terms in go_terms_per_entry:
        go_term_counts.update(terms)

    # Select top_n most frequent GO terms
    top_go_terms = [go for go, _ in go_term_counts.most_common(top_n)]

    # Initialize co-occurrence matrix with zeros
    cooc_matrix = pd.DataFrame(0, index=top_go_terms, columns=top_go_terms)

    # Count co-occurrences
    for terms in go_terms_per_entry:
        terms = set(terms).intersection(top_go_terms)  # Only consider top terms
        if include_diagonal:
            for term1, term2 in itertools.combinations_with_replacement(terms, 2):
                cooc_matrix.loc[term1, term2] += 1
                if term1 != term2:  # Don't double-count diagonal elements
                    cooc_matrix.loc[term2, term1] += 1
        else:
            for term1, term2 in itertools.combinations(terms, 2):
                cooc_matrix.loc[term1, term2] += 1
                cooc_matrix.loc[term2, term1] += 1

    # Rename indices and columns to include descriptions (if available)
    if go_term_descriptions:
        new_labels = {}
        for term in cooc_matrix.index:
            desc = go_term_descriptions.get(term, "")
            if desc:
                # Truncate long descriptions
                if len(desc) > 30:
                    desc = desc[:27] + "..."
                new_labels[term] = f"{term} ({desc})"
            else:
                new_labels[term] = term

        cooc_matrix.index = [new_labels[term] for term in cooc_matrix.index]
        cooc_matrix.columns = [new_labels[term] for term in cooc_matrix.columns]

    # Cluster terms if requested
    if cluster_terms:
        # Calculate distance matrix
        distance = hierarchy.distance.pdist(cooc_matrix.values)
        linkage = hierarchy.linkage(distance, method='average')

        # Compute the dendrogram
        dendro = hierarchy.dendrogram(linkage, no_plot=True)

        # Reorder the matrix based on dendrogram
        idx = dendro['leaves']
        cooc_matrix = cooc_matrix.iloc[idx, idx]

    # Plot heatmap
    plt.figure(figsize=figsize)

    # Determine colormap and scaling
    if log_scale and cooc_matrix.max().max() > 0:
        # Add a small value to avoid log(0)
        norm_data = np.log1p(cooc_matrix)
        cmap = "viridis"
    else:
        norm_data = cooc_matrix
        cmap = "Blues"

    # Determine whether to show annotations
    annot = include_annotations and (top_n <= 30)  # Only show annotations for smaller matrices

    # Create heatmap
    ax = sns.heatmap(
        norm_data,
        annot=annot,
        fmt=".0f" if not log_scale else ".1f",
        cmap=cmap,
        square=True,
        cbar_kws={"label": "Co-occurrence Count" if not log_scale else "Log(Co-occurrence Count + 1)"}
    )

    # Enhance plot appearance
    title_text = f"GO Term Co-occurrence Heatmap (Top {top_n} Terms)"
    if log_scale:
        title_text += " (Log Scale)"
    if include_diagonal:
        title_text += " with Self-occurrences"
    if cluster_terms:
        title_text += " (Clustered)"

    plt.title(title_text, fontsize=14)
    plt.xlabel("GO Terms")
    plt.ylabel("GO Terms")

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)

    # Add a grid to help with readability
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)

    plt.tight_layout()

    plt.show()


if __name__ == '__main__':

    json_file = "../01_raw_data/cleaned_raw_data/cleaned_uniprot_data.json"
    with open(json_file) as f:
        json_data = json.load(f)
        #GO_term_distribution_top_n(json_data, top_n=400)
        #GO_term_cooccurrence_heatmap(json_data, top_n=20)
        improved_GO_term_cooccurrence_heatmap(json_data, top_n=20, include_diagonal=False)
