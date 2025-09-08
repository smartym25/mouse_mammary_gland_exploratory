import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np 

filtered = pd.read_csv("data\GSE60450_filtered_metadata-1.csv")
genes = pd.read_csv("data\GSE60450_GeneLevel_NormalizedCPM.and_.TMM_data.csv")

genes.describe() #permette di vedere tutte le informazioni relative al dataframe 

filtered[["characteristics"]].iloc[1] #evidenzia la seconda riga del dataframe

#definiamo i nomi della prima colonna dei dataset 
df_filter = filtered.rename(columns={"Unnamed: 0":"sample_id"})

df_genes = genes.rename(columns={"Unnamed: 0":"gene_id"})

#creo la lista dei sample dei geni 
samples = [col for col in df_genes if col.startswith("GSM")]

seq_data = pd.melt(df_genes, id_vars=["gene_id", "gene_symbol"], value_vars=samples, var_name="sample", value_name="count")

#creo delle nuove colonne, le caratteristiche del gene, che si associano al suo id

fullinfo = pd.merge(seq_data, df_filter, left_on="sample", right_on="sample_id", how="outer")

allinfo = fullinfo.drop("sample_id", axis=1)

plt.figure(figsize=(18,12))

#visto che i count sono molto grandi di solito nei sequenziamenti viene utilizzato il log2 dei count 
#alcuni valore danno come origine ad -inf questo perchè il count era 0 di seguito aggiungiamo all'argomento dle logaritmo 1 
#questo non influenza in maniera eccessiva il risultato
allinfo["log2_count"] = np.log2(allinfo["count"]+1)

sns.set_theme(style="whitegrid")
sns.boxplot(allinfo, x="sample", y="log2_count") 

#eliminare i geni che non codificano la proteina o che non sono riconosciuti nel database swiss-prot
df_filtered_prot_gene = df[df["UniProtKB/Swiss-Prot ID"].notna()]

df_filtered_prot_gene.to_csv("data\prot_gene.csv", index=False)

df_filtered_prot_id = df_filtered_prot_gene["UniProtKB/Swiss-Prot ID"]

df_filtered_prot_id.to_csv("data\prot.csv", index=False, header=False)

#filtrare i geni togliendo quelli minori di tre, rimanendo solo quelli del prot.gene
#corrisponda tra sequenza peptidica dei fasta e le proteine nel datset prot_gene.csv


#togliere i rispettivi geni che non si evidenziano perchè il loro count è troppo basso log2 < 3 scartati

#eliminare le non proteine che codificano i geni 
