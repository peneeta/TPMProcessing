import pandas as pd

'''
Convert specified columns to integers, if none specified, converts the entire df to integers.
'''
def ConvertValuesToInts(df: pd.DataFrame, col_list: list[str] = None):
    
    if col_list:
        int_df = df[col_list].astype(int)
    else: int_df = df.astype(int)
    
    return int_df