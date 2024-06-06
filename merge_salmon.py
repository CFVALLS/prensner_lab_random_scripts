import pandas as pd
import os

current_directory = os.getcwd()
all_items = os.listdir(current_directory)

# Filter out only the directories
dirs = [item for item in all_items if os.path.isdir(os.path.join(current_directory, item))]
salmon_suffix = 'quant.sf'

df_merge = None
for directory in dirs:
    potential_path = os.path.join(current_directory, directory, salmon_suffix)
    if os.path.isfile(potential_path):
        file_path = potential_path
        
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        print(df.columns)
        
        # Rename the TPM column
        df.rename(columns={"TPM": directory}, inplace=True)
        print(df.head())
        
        # Merge DataFrames
        if df_merge is not None:
            df_merge = pd.merge(left=df_merge, right=df[['Name', directory]], how='outer', on='Name')
        else:
            df_merge = df[['Name', directory]]

        print(df_merge.shape)
if df_merge is not None:
    print(df_merge.head())
else:
    print("No quant.sf files found in any directories.")

df_merge.to_csv('salmon_merge.tsv', sep='\t', index=False)
print(df_merge.head())