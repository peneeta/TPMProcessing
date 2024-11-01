import pandas as pd

'''
SelectSamples filters a df based on user-defined samples in a string array format. Default: sample_list = None
df: a counts or tpm matrix where row names are genes and col names are samples
sample_list: an array of strings representing samples to select
'''
def SelectSamples(df: pd.DataFrame, sample_list: list[str] = None):
    #if !isinstance(tpm_df, pd.DataFrame):
    if sample_list:
        subset = df.iloc[:, sample_list]
        return subset
    
    else: return df

            
        
    