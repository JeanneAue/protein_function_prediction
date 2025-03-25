import json
import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (f1_score, precision_score, recall_score,
                           hamming_loss, accuracy_score,
                           average_precision_score, roc_auc_score,
                           precision_recall_curve, roc_curve,
                           classification_report, confusion_matrix,
                           multilabel_confusion_matrix)
import pandas as pd
import joblib
import xgboost
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MultiLabelBinarizer


def get_parents(term_id, hierarchy):
    parents_of_term = set()
    if hierarchy:
        for element in hierarchy:
            if element:
                if element["id"] == term_id:
                    for parent in element.get("parents", []):
                        if parent:
                            parent_id = parent.get("id")
                            parents_of_term.add(parent_id)
                else:
                    parent_ids = get_parents(term_id, element["parents"])
                    parents_of_term.update(parent_ids)
    return parents_of_term


def extract_go_parents_recursive(go_hierarchies, term_id, collected_parents=None):
    if collected_parents is None:
        collected_parents = set()
    collected_parents.add(term_id)
    parent_ids = get_parents(term_id, go_hierarchies)
    for parent_id in parent_ids:
        if parent_id and parent_id not in collected_parents:
            extract_go_parents_recursive(go_hierarchies, parent_id, collected_parents)
    return collected_parents


def retrieve_test_data():
    """Retrieve test data for evaluation."""
    # check if files x_test.joblib and y_test.joblib exist
    if os.path.exists('./x_test.joblib') and os.path.exists('./y_test.joblib'):
        print("Test data found on disk. Loading...")
        return joblib.load('./x_test.joblib'), joblib.load('./y_test.joblib')

    with open("../03_embeddings/protein_data_with_embeddings_and_hierarchy.json", "r") as infile:
        data = json.load(infile)

    X = np.array([entry["embedding"] for entry in data])
    y_raw = []

    for entry in data:
        go_terms = set(go["id"] for go in entry.get("go_annotations_with_ids", []))
        go_hierarchies = entry.get("go_hierarchies", [])
        all_terms = set()
        for term_id in go_terms:
            all_terms.update(extract_go_parents_recursive(go_hierarchies, term_id))
        y_raw.append(all_terms)

    mlb = MultiLabelBinarizer()
    y = mlb.fit_transform(y_raw)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    joblib.dump(X_test, './x_test.joblib')
    joblib.dump(y_test, './y_test.joblib')
    print("Test data saved to disk.")

    return X_test, y_test

# save x_test and y_test to a file
X_test, y_test = retrieve_test_data()

# Load your models and test data
# Replace these with your actual file paths
# load models from .pkl files as joblib # TODO Sarah
rf_filepath = "models/rf.pkl"
svm_filepath = "models/svm.pkl"
gb_filepath = "models/gb.pkl"
mlb_filepath = "models/mlb.pkl"

rf_model = joblib.load(rf_filepath)
svm_model = joblib.load(svm_filepath)
gb_model = joblib.load(gb_filepath)
mlb = joblib.load(mlb_filepath)  # If you have one

print("loaded models")

print("X_test shape:", X_test.shape)
print("y_test shape:", y_test.shape)



# Dictionary of models
models = {
    'Random Forest': rf_model,
    'SVM': svm_model,
    'Gradient Boosting': gb_model
}

# Color map for consistent colors across plots
colors = ['#f04729', '#f07c29', '#f0b429']
model_colors = {model_name: color for model_name, color in zip(models.keys(), colors)}

# Dictionary to store results
results = {}
y_preds = {}
y_probas = {}

# Get predictions from each model
for model_name, model in models.items():
    y_pred = model.predict(X_test)
    y_preds[model_name] = y_pred

    # Store probability predictions if available
    if hasattr(model, "predict_proba"):
        y_probas[model_name] = model.predict_proba(X_test)
    elif hasattr(model, "decision_function"):
        y_probas[model_name] = model.decision_function(X_test)
    else:
        y_probas[model_name] = None

    # Transform binary predictions to binary matrix
    y_preds_transformed = {}
    for model_name, y_pred in y_preds.items():
        # If y_pred is 1D, convert to binary matrix
        if y_pred.ndim == 1:
            y_preds_transformed[model_name] = mlb.transform(y_pred.reshape(-1, 1))
        else:
            y_preds_transformed[model_name] = y_pred

    # Calculate metrics
    results[model_name] = {
        'hamming_loss': hamming_loss(y_test, y_pred), # todo why did this fail with svm due to inconsistent shapes?
        'subset_accuracy': accuracy_score(y_test, y_pred, normalize=True),
        'micro_precision': precision_score(y_test, y_pred, average='micro'),
        'micro_recall': recall_score(y_test, y_pred, average='micro'),
        'micro_f1': f1_score(y_test, y_pred, average='micro'),
        'macro_precision': precision_score(y_test, y_pred, average='macro'),
        'macro_recall': recall_score(y_test, y_pred, average='macro'),
        'macro_f1': f1_score(y_test, y_pred, average='macro'),
        'per_class_f1': f1_score(y_test, y_pred, average=None),
        'per_class_precision': precision_score(y_test, y_pred, average=None),
        'per_class_recall': recall_score(y_test, y_pred, average=None)
    }
    print("finished model:", model_name)
    print("results:", results[model_name])

# Convert results to DataFrame for easier plotting
df_results = pd.DataFrame(results).T
print("gathered results. continuing...")

# 1. Overall Model Comparison Plot
def plot_overall_comparison(df_results):
    """Plot overall model comparison metrics."""
    metrics = ['micro_f1', 'macro_f1', 'micro_precision', 'micro_recall',
               'macro_precision', 'macro_recall']

    plt.figure(figsize=(12, 6))

    # Create a grouped bar chart
    x_pos = np.arange(len(metrics))
    width = 0.25
    offsets = [-width, 0, width]

    for i, (model_name, color) in enumerate(model_colors.items()):
        values = [df_results.loc[model_name, metric] for metric in metrics]
        plt.bar(x_pos + offsets[i], values, width, label=model_name, color=color, alpha=0.8)

    plt.xlabel('Metrics')
    plt.ylabel('Score')
    plt.title('Model Performance Comparison')
    plt.xticks(x_pos, metrics)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

# 2. Hamming Loss Plot (lower is better)
def plot_hamming_loss(df_results):
    """Plot hamming loss for each model (lower is better)."""
    plt.figure(figsize=(8, 5))
    plt.bar(df_results.index, df_results['hamming_loss'],
            color=[model_colors[model] for model in df_results.index], alpha=0.8)
    plt.xlabel('Model')
    plt.ylabel('Hamming Loss')
    plt.title('Hamming Loss (lower is better)')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.ylim(0, min(1.0, df_results['hamming_loss'].max() * 1.5))
    plt.tight_layout()
    plt.show()

# 3. Per-class F1 Score Plot
def plot_per_class_f1(results, mlb):
    """Plot F1 score for each class and each model."""
    # Get class names
    class_names = mlb.classes_ if mlb is not None else [f"Class {i}" for i in range(len(results['Random Forest']['per_class_f1']))]

    # Create a DataFrame with per-class F1 scores
    data = {}
    for model_name in results:
        data[model_name] = results[model_name]['per_class_f1']

    df_f1 = pd.DataFrame(data, index=class_names)

    # Select top N classes based on average F1 score
    n_classes = min(20, len(df_f1))  # Limit to 20 classes for readability
    top_classes = df_f1.mean(axis=1).sort_values(ascending=False).index[:n_classes]
    df_f1_top = df_f1.loc[top_classes]

    plt.figure(figsize=(14, 8))
    df_f1_top.plot(kind='bar', rot=90, figsize=(14, 8))
    plt.title(f'Per-class F1 Score (Top {n_classes} classes)')
    plt.xlabel('GO Term')
    plt.ylabel('F1 Score')
    plt.legend(title='Model')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

    return df_f1

# 4. Micro-Averaged ROC Curve
def plot_micro_roc(y_test, y_probas, model_colors):
    """Plot micro-averaged ROC curve for each model."""
    plt.figure(figsize=(10, 8))

    for model_name, y_proba in y_probas.items():
        if y_proba is not None:  # Skip models without probability predictions
            # Compute ROC curve and ROC area for each class
            fpr = dict()
            tpr = dict()
            roc_auc = dict()

            n_classes = y_test.shape[1]

            # Compute micro-average ROC curve and ROC area
            fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_proba.ravel())
            roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

            plt.plot(fpr["micro"], tpr["micro"],
                     label=f'{model_name} (AUC = {roc_auc["micro"]:.3f})',
                     color=model_colors[model_name], lw=2)

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Micro-Averaged ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

# 5. Precision-Recall Curve
def plot_precision_recall(y_test, y_probas, model_colors):
    """Plot precision-recall curve for each model."""
    plt.figure(figsize=(10, 8))

    for model_name, y_proba in y_probas.items():
        if y_proba is not None:  # Skip models without probability predictions
            # Compute average precision
            avg_precision = average_precision_score(y_test.ravel(), y_proba.ravel())

            # Compute precision-recall curve
            precision, recall, _ = precision_recall_curve(y_test.ravel(), y_proba.ravel())

            plt.plot(recall, precision,
                     label=f'{model_name} (AP = {avg_precision:.3f})',
                     color=model_colors[model_name], lw=2)

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc="best")
    plt.grid(linestyle='--', alpha=0.7)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.tight_layout()
    plt.show()

# 6. Class frequency vs. F1 score plot
def plot_class_freq_vs_f1(y_test, results, mlb):
    """Plot the relationship between class frequency and F1 score."""
    # Get class frequencies
    class_freq = y_test.sum(axis=0)
    class_names = mlb.classes_ if mlb is not None else [f"Class {i}" for i in range(len(class_freq))]

    plt.figure(figsize=(12, 8))

    for model_name, color in model_colors.items():
        f1_scores = results[model_name]['per_class_f1']
        plt.scatter(class_freq, f1_scores, alpha=0.6, label=model_name, color=color)

    plt.xscale('log')
    plt.xlabel('Class Frequency (log scale)')
    plt.ylabel('F1 Score')
    plt.title('Class Frequency vs. F1 Score')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# 7. Confusion matrix heatmap for top N classes
def plot_top_classes_cm(y_test, y_pred, mlb, n=10):
    """Plot confusion matrices for top N classes."""
    # Get class names
    class_names = mlb.classes_ if mlb is not None else [f"Class {i}" for i in range(y_test.shape[1])]

    # Get top N classes by frequency
    class_freq = y_test.sum(axis=0)
    top_indices = np.argsort(-class_freq)[:n]
    top_class_names = [class_names[i] for i in top_indices]

    # Get confusion matrices for these classes
    cm = multilabel_confusion_matrix(y_test, y_pred)

    # Plot confusion matrices for top classes
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()

    for i, idx in enumerate(top_indices[:n]):
        sns.heatmap(cm[idx], annot=True, fmt="d", cmap="Blues", ax=axes[i])
        axes[i].set_title(f"{top_class_names[i]}")
        axes[i].set_xlabel('')
        axes[i].set_ylabel('')

    plt.tight_layout()
    plt.show()

# Example usage
plot_overall_comparison(df_results)
plot_hamming_loss(df_results)
plot_per_class_f1(results, mlb)

# For models with probability outputs:
if all(prob is not None for prob in y_probas.values()):
    plot_micro_roc(y_test, y_probas, model_colors)
    plot_precision_recall(y_test, y_probas, model_colors)

plot_class_freq_vs_f1(y_test, results, mlb)

# For a selected model (e.g., best performing one)
best_model = max(results, key=lambda x: results[x]['micro_f1'])
plot_top_classes_cm(y_test, y_preds[best_model], mlb)

# Create a results summary table
print("\nModel Performance Summary:")
print(df_results[['micro_f1', 'macro_f1', 'hamming_loss', 'subset_accuracy']].round(3))

# Export results to CSV
df_results.to_csv('model_comparison_results.csv')