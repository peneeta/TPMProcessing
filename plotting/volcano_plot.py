import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import adjust_text

# Plots a volcano plot using deseq_df
def PlotVolcano(deseq_df: pd.DataFrame, outputFile: str, plot_title: str, lbl_logFDR_threshold: float, lbl_FC_threshold: float):
    # Create a scatter plot of p-adjusted vs log2fc
    plt.scatter(x=deseq_df['log2FoldChange'], y=deseq_df['padj'].apply(lambda x:-np.log10(x)), s=1)

    # Get down and upregulated genes
    down, up = GetDownAndUpregulatedGenes(deseq_df)

    # Scatter up and downregulated genes
    plt.scatter(x=down['log2FoldChange'],y=down['padj'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    plt.scatter(x=up['log2FoldChange'],y=up['padj'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

    # Add labels to most significant genes (both down and upregulated)
    lbls = GetLabels(down, up, lbl_logFDR_threshold, lbl_FC_threshold)
    
    # Put labels on figure
    adjust_text(lbls,
                expand=(2, 3),
                arrowprops=dict(arrowstyle='-', connectionstyle='arc3,rad=0', lw=0.5))

    # Labels for plot
    plt.xlabel("logFC")
    plt.ylabel("-logFDR")
    plt.title(plot_title)
    plt.axvline(-2,color="grey",linestyle="--")
    plt.axvline(2,color="grey",linestyle="--")
    plt.axhline(2,color="grey",linestyle="--")
    plt.legend()
    plt.tight_layout()

    # Save plot as jpg
    plt.savefig(os.path.join(outputFile+'.jpg'))


def GetDownAndUpregulatedGenes(deseq_df: pd.DataFrame):
    down = deseq_df[(deseq_df['log2FoldChange']<=-2)&(deseq_df['padj']<=0.01)&(deseq_df['padj']!=0)] # remove p-values = 0
    
    up = deseq_df[(deseq_df['log2FoldChange']>=2)&(deseq_df['padj']<=0.01)&(deseq_df['padj']!=0)]
    
    return down, up

    
def GetLabels(down_df, up_df, log_threshold, fc_threshold):
    lbls=[]
    
    # Get downregulated labels
    down_labels = down_df[(down_df['log2FoldChange']<=-log_threshold)|((down_df['padj'].apply(lambda x:-np.log10(x)))>fc_threshold)]
    
    # Get upregulated labels
    up_labels = up_df[(up_df['log2FoldChange']>=log_threshold)|((up['padj'].apply(lambda x:-np.log10(x)))>fc_threshold)]
    
    for i,r in up_labels.iterrows():
        lbls.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=r['gene_name']))
    
    
    for i,r in down_labels.iterrows():
        lbls.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=r['gene_name']))
    
    return lbls