# miRNA-embedding-project
MatureProfilesConstruction - Converts raw profiles to profiles with relevant data.
DocumentsReduction - Responsible of the conversion process from profiles (after recreation in MatureProfilesConstruction) to text documents. 
WordEmbedding - Reads the text documents and creates the embeddings for the miRNAs. 
Classification - Checking the quiality of the embedding by classifying the profiles to normal and tumor, before and after the embedding process. 
TargetGenes - Responsible of computing the similarity between miRNAs according to their overlap in target genes.
MiRNADiseases - Responsible of computing the similarity between miRNAs according to their overlap in associated diseases.
MiRNACorrelation - Responsible of computing the similarity between miRNAs according to the correlation between their counts in the profiles.
MergeDatasets - Helps merging the datasets with the different similarity metrics.
CorrelationBetweenMetrics - A Google Colab notebook with measurements of the correlation between the similarity metrics. 
