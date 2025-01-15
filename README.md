# protein_function_prediction

The raw_data folder should contain 3 files:
- uniprotkb_human_AND_model_organism_9606_2025_01_14.fasta
- uniprotkb_human_AND_model_organism_9606_2025_01_14.xml
- ...

Download available at: https://www.uniprot.org/uniprotkb?query=human&facets=model_organism:9606,reviewed:true


GO-term quick view: https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=P01308

## Fragen:
- embeddings/padding -> "PAD" padding in Ordung oder manuell über BED files erstellen (wenn ja woher?)
  - Antwort: PAD reicht!!
  - TODO: model.eval ?
  - TODO: batching
- -welche Infos pro Protein sollen wir sinnvollerweise speichern? (aktuell F:...)-
- welche Plots/Infos sollen wir bei der Analyse den Fokus drauf legen?
  - verteilung mit allen hierarchien
  - das war die Verteilung, so haben wir entschieden, wo wir den cutoff machen
  - wie groß sind embeddings (jeweils für batching!!)
- vorhersage auch für hierarchie der GO terms oder nur konkrete GO terms
  - in den Rohdaten haben wir mehrere Funktionen
- embeddings: nicht zu viel Zeit drauf verwenden (es gibt nicht so viele Cluster, weil zu hochdimensional)
  - auf jeden Fall mal plotten
  - andere Formen für Analyse