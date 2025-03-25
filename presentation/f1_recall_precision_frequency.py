import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# HPI Colors
HPI_COLORS = {
    "Random Forest": "#E3001B",  # HPI Red
    "SVM": "#F37521",            # HPI Orange
    "Gradient Boosting": "#FFCC00"  # HPI Yellow
}

# Helper function to extract classification report metrics from log files
def extract_classification_report(log_text):
    # Regular expression for your specific format
    pattern = r"(GO:\d+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)"
    matches = re.findall(pattern, log_text)

    if not matches:
        print("No matches found. Please check the log file format.")
        return pd.DataFrame()

    # Convert matches to a DataFrame
    report_data = []
    for match in matches:
        label, precision, recall, f1, support = match
        report_data.append({
            "Label": label,
            "Precision": float(precision),
            "Recall": float(recall),
            "F1-Score": float(f1),
            "Support": int(support)
        })

    # Filter out rows with zero support
    df = pd.DataFrame(report_data)
    df = df[df["Support"] > 0]
    return df

# Function to group GO terms by frequency bins
def group_by_frequency(report, bins):
    # Add a "Frequency Bin" column based on support bins
    report["Frequency Bin"] = pd.cut(report["Support"], bins=bins, labels=bins[:-1])
    
    # Group by "Frequency Bin" and calculate mean only for numeric columns
    numeric_columns = ["Precision", "Recall", "F1-Score"]
    grouped = report.groupby("Frequency Bin")[numeric_columns].mean()
    
    return grouped

# Load logs (replace paths with actual file paths)
log_files = {
    "Random Forest": "./reports/rf.log",
    "SVM": "./reports/bio-medical-21074577.log",
    "Gradient Boosting": "./reports/bio-medical-31074578.log"
}

classification_reports = {}

# Extract classification reports from each log file
for model_name, log_path in log_files.items():
    with open(log_path, 'r') as file:
        log_content = file.read()
        classification_reports[model_name] = extract_classification_report(log_content)

# Define bins for frequency grouping using log spacing for better granularity
frequency_bins = np.logspace(0, np.log10(800), num=10).astype(int)
grouped_metrics = {}

for model_name, report in classification_reports.items():
    grouped_metrics[model_name] = group_by_frequency(report, frequency_bins)

# Plot F1-Score by frequency bins using HPI colors
plt.figure(figsize=(14, 7))
for model_name, metrics in grouped_metrics.items():
    plt.plot(metrics.index, metrics["F1-Score"], label=f"F1-Score ({model_name})", marker="o", color=HPI_COLORS[model_name])
plt.title("F1-Score by GO Term Frequency", fontsize=14)
plt.xlabel("GO Term Frequency (Support Bins)", fontsize=12)
plt.ylabel("Average F1-Score", fontsize=12)
plt.legend()
plt.grid()
plt.show()

# Plot Precision by frequency bins using HPI colors
plt.figure(figsize=(14, 7))
for model_name, metrics in grouped_metrics.items():
    plt.plot(metrics.index, metrics["Precision"], label=f"Precision ({model_name})", marker="o", color=HPI_COLORS[model_name])
plt.title("Precision by GO Term Frequency", fontsize=14)
plt.xlabel("GO Term Frequency (Support Bins)", fontsize=12)
plt.ylabel("Average Precision", fontsize=12)
plt.legend()
plt.grid()
plt.show()

# Plot Recall by frequency bins using HPI colors
plt.figure(figsize=(14, 7))
for model_name, metrics in grouped_metrics.items():
    plt.plot(metrics.index, metrics["Recall"], label=f"Recall ({model_name})", marker="o", color=HPI_COLORS[model_name])
plt.title("Recall by GO Term Frequency", fontsize=14)
plt.xlabel("GO Term Frequency (Support Bins)", fontsize=12)
plt.ylabel("Average Recall", fontsize=12)
plt.legend()
plt.grid()
plt.show()
