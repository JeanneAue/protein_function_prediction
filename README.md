# Protein Function Prediction from sequence data

Sources for exploration:

UniProt - Download available at: https://www.uniprot.org/uniprotkb?query=human&facets=model_organism:9606,reviewed:true

QuickGo - GO-term quick lookup: https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=P01308

GOA - Gene Ontology Annotation: https://www.ebi.ac.uk/GOA/

## Prepare data for cluster 
1. **Data Collection**: Download human proteins from UniProt **or**: unzip `01_raw_data/uniprotkb_human_AND_model_organism_9606_2025_01_14.json.zip`
2. **Clean Data**: 

   a) execute `01_raw_data/clean_data.py` and 

   b) execute `01_raw_data/parsing_for_embedding.py`
3. **Add GO terms**: Execute `02_data_analysis/filter_and_add_GO_terms_json.py`
3. **Create embeddings using BERT**: Execute `03_embeddings/embedding_cluster.py`
4. **Add GO term hierarchie**: Execute `02_data_analysis/go_term_add_hierarchie.py`

## Cluster execution
1. Copy `03_embeddings/protein_data_with_embeddings_and_hierarchy.json` and one (or multiple) of the three models to the cluster node 

   a) `04_prediction/gradient_boosting_cluster.py`

   b) `04_prediction/random_forest_cluster.py`

   c) `04_prediction/support_vector_machine_cluster.py`
2. Execute the model on the cluster node

Troubleshooting: Look at the logs in reports folder for library versions and configurations.

## Inference
Try out the notebooks in `05_inference`
