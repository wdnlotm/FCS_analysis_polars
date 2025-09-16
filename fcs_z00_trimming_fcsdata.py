"""
python fcs_z00_trimming_fcsdata.py --dataloc $FCS_DATA_DIR
or in Jupyter Notebook
%run fcs_z00_trimming_fcsdata.py --dataloc $FCS_DATA_DIR
"""
import os, argparse, yaml
import gzip
import flowkit as fk
import numpy as np
import pandas as pd
import polars as pl

# function to add sample id to fcs file data
def augmented_df(df, sid):
    # df: polars dataframe
    return df.with_columns(pl.lit(sid).alias(  "('sample_id', 'sample_id')"  ))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FCS trimming. take datadir')
    parser.add_argument('-dl','--dataloc', type=str, default="./", help='data location')
    args, unknown = parser.parse_known_args()
    
    print(os.getcwd())
    
    path_files = args.dataloc #"data_cytof/"
    
    samples_all = fk.load_samples(path_files, filename_as_id=True)
    print(f'There are {len(samples_all)} fcs files.')
    n = 7
    print(f'First {n}: {samples_all[0:n]}')
    
    sample_df_list = [augmented_df( pl.from_pandas(samples_all[ii].as_dataframe(source='raw')), 
                                    samples_all[ii].id) for ii in range(len(samples_all))]
    fcs_df=pl.concat(sample_df_list)
    print(fcs_df.shape, fcs_df.shape[0]/1e6)
    
    # This is for polars dataframe
    col_orig=fcs_df.columns
    col_orig_list=list(col_orig)
    
    col_name_list2=[]
    # marker_col=[]
    
    for ii in range(len(col_orig_list)):
        loc1 = col_orig_list[ii].index('(') #remove hypen
        loc2 = col_orig_list[ii].index(',')
        loc3 = col_orig_list[ii].index(')')
    
        col_name_list2.append( (  col_orig_list[ii][(loc1+2):(loc2-1)]  , 
                                  col_orig_list[ii][(loc2+3):(loc3-1)].replace("-","")   )  )
    
    new_columns=pd.MultiIndex.from_arrays([[x[0] for x in col_name_list2],[y[1] for y in col_name_list2]], 
                              names=['pnn', 'pns'])
    
    new_columns_single = [x[1] for x in new_columns]
    # new_columns_single
    fcs_df.columns=new_columns_single #single column names
    
    new_columns_single_df = pd.DataFrame({'antigen': new_columns_single})
    new_columns_single_df.to_csv('new_columns.csv', index=False)    #save this csv and use it to mark which columns to include
                                                       # return "columns_2_use.csv" 
    
    print(f'\nCheck out new_columns.csv file. Add a column named "use" to indicate which one to use (1) and which one not (0).')
    print(f'Save the file as columns_2_use.csv')
    print('\n')
    
    restart = input('Type anything here when columns_2_use_csv is ready.')
    
    columns_2_use_df = pd.read_csv('columns_2_use.csv')
    
    fcs_df_trimmed = fcs_df[list(columns_2_use_df.loc[columns_2_use_df['use']==1, 'antigen'])]
    
    # further simplify column names
    column_list = list(fcs_df_trimmed.columns)
    new_column_list = [cn[(cn.index('_')+1):(len(cn))].replace('_',"") for cn in column_list[0:(fcs_df_trimmed.shape[1]-1)]]
    fcs_df_trimmed.columns = new_column_list+['sample_id']
    fcs_df_trimmed.head(2)
    
    file_path = "fcs_df_trimmed.csv.gz"
    with gzip.open(file_path, "wb") as f:
        fcs_df_trimmed.write_csv(f)
    
    print(f'\n\nfcs data only with relavant markers saved. {file_path}')

