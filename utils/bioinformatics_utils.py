from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import pandas as pd

'''
ConductDeseq uses pydeseq2 to conduct deseq2 using a user-specified counts matrix, metadata, and design factor.
count_mtrx: a pandas dataframe 
metadata_df: a dataframe matching each sample to a type (i.e. if a sample is named WT_Replicate2, its cell type metadata would be 'WT' or 'wildtype')
design_fct: a string representing the comparison Deseq2 should make. This should match the column or columns in metadata_df
deseq_contrast: an array of strings representing the contrasts of interest to be passed into the DeseqStats function
logFC_threshold: default = 2, return genes that are greater than this threshold
pVal_threshold: default = 0.01, return genes that are greater than this threshold
'''
def ConductDeseq(count_mtrx: pd.DataFrame, metadata_df: pd.DataFrame, design_fct: str, deseq_contrast: list[str], logFC_threshold: float = 2, pVal_threshold: float = 0.01):
    # Transpose counts so samples are the row values (need this format for deseq)
    count_mtrx_t = count_mtrx.T
    
    dds = DeseqDataSet(counts=count_mtrx_t, metadata=metadata_df, design_factors = design_fct)
    dds.deseq2()

    # Distinguish between KO vs WT samples 
    stat_res = DeseqStats(dds, contrast = (design_fct, deseq_contrast))
    
    stat_res.summary() # must call this to access results_df
    res = stat_res.results_df

    # Reset index and set gene names as regular column
    res.reset_index(inplace=True)
    res = res.rename(columns={"index": "gene_name"})

    # DE genes that are significant
    sigs = res[(res.padj < pVal_threshold)& (abs(res.log2FoldChange)>logFC_threshold)]
    
    # return
    return res, sigs

