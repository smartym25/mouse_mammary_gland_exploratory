import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np 
from Bio import SeqIO

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

#FILTRO LA TABELLA DEI GENI IN BASE AD UN SAMPLE
#selezione il sample
df_sample2 = allinfo[allinfo["sample"] == "GSM1480292"]

#tolgo i geni che hanno il log2 < 2
df_sample2_filt = df_sample2[df_sample2["log2_count"] > 2]

#ordiniamo in ordine decrescente log2, il valore più alto indica il gene più espressivo
df_order_log2 = df_sample2_filt.sort_values(by="log2_count", ascending=False)

#salvo il file
df_order_log2.to_csv("data\sample2.csv", index=False)

#dopo aver filtrato la tabella, creiamo un file con solo le proteine da mettere in ensembl
df_sample2_gene_id = df_order_log2["gene_id"]

df_sample2_gene_id.to_csv("data\prot.csv", index=False, header=False)

#otteniamo il dataset da ensembl e lo nominiamo prot_genes.csv
df_prot_genes = pd.read_csv("data\prot_genes.csv")

df_prot_genes_rename = df_prot_genes.rename(columns={"Gene stable ID":"gene_id",
                                                     "Gene name":"gene_symbol"})

#salvo il file con i geni che codificano proteine 
df_prot_genes_rename.to_csv("data\prot_gene_filter.csv", index=False)

#unire prot_gene_filter.csv e sample2.csv
df_prot_gene_final = pd.merge(df_order_log2, 
    df_prot_genes_rename[["gene_id", "UniProtKB/Swiss-Prot ID"]],
    on="gene_id", 
    how="left")

#eliminare i geni che non codificano la proteina o che non sono riconosciuti nel database swiss-prot
df = df_prot_gene_final[df_prot_gene_final["UniProtKB/Swiss-Prot ID"].notnull()]

df_final = df.drop_duplicates()

df_final.to_csv("data\prot_genes_final.csv", index=False)

#ATTRIBUIRE AD UN SWISS ID LA SEQUENZA PEPTIDICA DELLA PROTEINA 
#ricavo in una lista i swiss id delle proteine 
df_prot_id = df_final["UniProtKB/Swiss-Prot ID"]

df_prot_id.to_csv("data\prot.csv", index=False, header=False)

#utilizziamo Biopython per leggere il file fasta contenente la seuquenza peptidica
def transform_list(filename, option="id"):
    
    results = []
    
    if option not in ("id", "seq"):
        raise ValueError("option must be id or seq")
    
    for record in SeqIO.parse(filename, "fasta"):
        if option == "id":
            results.append(record.id[:18])
        elif option == "seq":
            results.append(record.seq)
    return results

df_seq_prot = pd.DataFrame({"gene_id": transform_list("fasta_prot.fasta"),
                            "seq_prot": transform_list("fasta_prot.fasta", "seq")})

#uniamo la tabella con le sequenze peptidiche e df_final 
table = pd.merge(df_final, 
    df_seq_prot[["gene_id", "seq_prot"]],
    on="gene_id", 
    how="left")

table_2 = table.drop_duplicates()

table_3 = table_2.drop(["characteristics", "immunophenotype", "developmental stage", "sample"], axis=1)

table_3.to_csv("final_info.csv", index=False)
