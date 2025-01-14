# Protein Function Prediction from Protein Sequence

## 1. Data Acquisition
- **Download protein sequences and annotations**:
  - Download all **human proteins** from **UniProt**.
  - Get **Gene Ontology (GO) annotations** for these proteins, classified into (for us most likely only molecular function is relevant):
    - **Molecular Function**
    - **Cellular Component**
    - **Biological Process**


## 2. Data Preparation
- Perform a **descriptive analysis** of the downloaded data, including:
  - Distribution of **sequence lengths**.
  - Distribution of **GO annotations**.
  - Perform **data cleaning**:
    - Handle missing annotations.
    - Remove redundant data.
- Decide if further preprocessing of protein sequences is required, such as tokenization for embedding models.


## 3. Protein Embedding
- Use a **protein language model** to embed the protein sequences into numerical representations. Examples include:
  - **ProtTrans** models like `ProtBERT` or `ProtT5`.
  - **ESM (Evolutionary Scale Modeling)**.
  - **SeqVec** (based on ELMo).
- The embeddings will serve as input features for machine learning.


## 4. Train a Classifier
- Use **scikit-learn** to build a classification model:
  - **Input**: Protein embeddings (features).
  - **Output**: Predicted GO annotations (target labels).
- Example classifiers:
  - Logistic Regression.
  - Random Forest.
  - Support Vector Machines (SVM).
  - Gradient Boosting (e.g., XGBoost or LightGBM).



## 5. Evaluation
- Compare the classifier's predictions with the actual GO annotations to assess performance using metrics such as:
  - **Accuracy**
  - **Precision**, **Recall**, **F1-Score**
  - **ROC-AUC** (if applicable).
- Compare your results with **other approaches** from the literature.



## 6. Benchmarking
- Review existing literature on **protein function prediction** and compare your approach with:
  - Reported metrics.
  - Models used (e.g., deep learning vs traditional ML).
  - Embedding techniques.


# Data cleaning

How is our data structured

* Basically everything is inside of a `results` array. This array ontains objects which are the entries each starting with entryType.
* For each entry there is a `uniProtKBCrossReference` key we want to keep an entry if inside in one of the `properties` fields there is an object with

```json
"properties": [
{
  "key": "GoTerm",
  "value": "F:xxxxx"
},{
}
]
```

* What we want to remove on the level of `entryType` and `uniProtKBCrossReference`
  * entryAudit
  * organism
  * references
