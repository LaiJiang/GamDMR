import pandas as pd

# Define input and output file paths
input_file = '/mnt/c/Per/LaiJiang/Project/UQAC/meth/dat/peek_methylation_sva.csv'
#output_file = '/mnt/c/Per/LaiJiang/Project/UQAC/meth/dat/peek_methylation_sva_short.csv' with 0.05
output_file = '/mnt/c/Per/LaiJiang/Project/UQAC/meth/dat/peek_methylation_sva_0035.csv' #with 0.035

# Set an appropriate chunksize (number of rows per chunk)
chunksize = 10000

# This flag helps us write the header only for the first chunk
first_chunk = True

# Process the file chunk by chunk
for chunk in pd.read_csv(input_file, sep="\t", chunksize=chunksize):
    # Identify methylation columns (those that end with '_meth')
    meth_cols = [col for col in chunk.columns if col.endswith('_meth')]
    
    # Calculate the variance across the methylation columns for each row
    # (skipna=True handles any missing values)
    chunk['var_meth'] = chunk[meth_cols].var(axis=1, skipna=True)
    
    # Filter rows where variance 
    filtered_chunk = chunk[chunk['var_meth'] >= 0.035]
    
    # Drop the temporary variance column before writing out
    filtered_chunk = filtered_chunk.drop(columns=['var_meth'])
    
    # Write the filtered rows to the output file
    if first_chunk:
        filtered_chunk.to_csv(output_file, sep="\t", index=False, mode='w')
        first_chunk = False
    else:
        filtered_chunk.to_csv(output_file, sep="\t", index=False, mode='a', header=False)

print("Processing complete. Filtered data saved to:", output_file)
