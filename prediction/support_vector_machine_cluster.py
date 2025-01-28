import json
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.multioutput import MultiOutputClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report, hamming_loss
import joblib

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

def save_checkpoint(model, mlb, filename="checkpoint.pkl"):
    print(f"Saving checkpoint to {filename}...")
    joblib.dump({"model": model, "mlb": mlb}, filename)
    print("Checkpoint saved.")

def load_checkpoint(filename="checkpoint.pkl"):
    if os.path.exists(filename):
        print(f"Loading checkpoint from {filename}...")
        checkpoint = joblib.load(filename)
        print("Checkpoint loaded.")
        return checkpoint["model"], checkpoint["mlb"]
    return None, None


with open("protein_data_with_embeddings_and_hierarchy.json", "r") as infile:
    data = json.load(infile)

X = np.array([entry["embedding"] for entry in data])
y_raw = []


print("Processing annotations and hierarchies...")
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

# Identify and filter out columns with only one class in the training set
non_constant_columns = []
for i in range(y_train.shape[1]):
    unique_vals = np.unique(y_train[:, i])
    if len(unique_vals) > 1:  # means there's both 0 and 1
        non_constant_columns.append(i)

non_constant_columns = np.array(non_constant_columns)
if len(non_constant_columns) == 0:
    raise ValueError("All labels have only one class in the training set.")

y_train_filtered = y_train[:, non_constant_columns]
y_test_filtered = y_test[:, non_constant_columns]

classifier, mlb_checkpoint = load_checkpoint()
if classifier is None:
    print("No checkpoint found. Starting fresh...")
    base_model = SVC(probability=True, random_state=42)
    classifier = MultiOutputClassifier(base_model)
else:
    print("Resuming training from checkpoint...")

try:
    print("Training the classifier...")
    classifier.fit(X_train, y_train_filtered)
    print("Training completed.")

    # Save checkpoint
    save_checkpoint(classifier, mlb, filename="svm_protein_function_checkpoint.pkl")

    print("Evaluating the classifier...")
    y_pred_filtered = classifier.predict(X_test)
    used_label_names = [mlb.classes_[i] for i in non_constant_columns]

    print("Classification Report (for non-constant labels only):")
    print(classification_report(y_test_filtered, y_pred_filtered, target_names=used_label_names))

    print("Hamming Loss:", hamming_loss(y_test_filtered, y_pred_filtered))

    joblib.dump(classifier, "svm_protein_function_classifier_with_hierarchy.pkl")
    joblib.dump(mlb, "svm_go_mlb_with_hierarchy.pkl")
    print("Model and label binarizer saved.")

except KeyboardInterrupt:
    print("Training interrupted. Saving checkpoint...")
    save_checkpoint(classifier, mlb, filename="svm_protein_function_checkpoint.pkl")
except Exception as e:
    print(f"An error occurred: {e}. Saving checkpoint...")
    save_checkpoint(classifier, mlb, filename="svm_protein_function_checkpoint.pkl")
