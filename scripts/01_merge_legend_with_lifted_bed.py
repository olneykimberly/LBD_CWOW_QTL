import pandas as pd

# Read the output.bed file
bed_df = pd.read_csv('1000GP_Phase3_combined_lifted_to_GRCh38.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'id'])

# Read the 1000GP_Phase3_combined.legend file with appropriate header and data types
legend_dtypes = {
    'id': str,
    'chr': str,
    'position': int,
    'a0': str,
    'a1': str,
    'TYPE': str,
    'AFR': float,
    'AMR': float,
    'EAS': float,
    'EUR': float,
    'SAS': float,
    'ALL': float
}

legend_df = pd.read_csv('1000GP_Phase3_combined.legend', sep=' ', dtype=legend_dtypes, header=0)

# Merge the dataframes on the 'id' column from both files
merged_df = pd.merge(legend_df, bed_df, on='id', how='left')

# Save the result to a new file
merged_df.to_csv('1000GP_Phase3_combined_lifted_to_GRCh38_merged_output.txt', sep='\t', index=False)

print("Merging completed. The output file is '1000GP_Phase3_combined_lifted_to_GRCh38_merged_output.txt'.")

