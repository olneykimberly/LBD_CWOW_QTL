import pandas as pd

# Specify data types for each column. Set problematic columns to float.
dtype_dict = {
    'id': str,
    'chr_x': str,
    'position': str,
    'a0': str,
    'a1': str,
    'TYPE': str,
    'AFR': float,
    'AMR': float,
    'EAS': float,
    'EUR': float,
    'SAS': float,
    'ALL': float,
    'chr_y': str,
    'start': float,  # Change this to float64 to handle mixed types
    'end': float     # This should also be float64 to handle any decimal points
}

# Read the merged output file with specified data types and low_memory=False
merged_df = pd.read_csv('1000GP_Phase3_combined_lifted_to_GRCh38_merged_output.txt', sep='\t', dtype=dtype_dict, low_memory=False)

# Rename columns
merged_df.rename(columns={'chr_y': 'chr'}, inplace=True)

# Convert the 'end' column to string without .0
merged_df['end'] = merged_df['end'].astype(str).str.replace('.0', '', regex=False)

# Split the 'id' column by ":" and replace the value between the first and second ":" with the 'end' value
def reformat_id(row):
    parts = row['id'].split(':')
    parts[1] = row['end']
    return ':'.join(parts)

merged_df['id'] = merged_df.apply(reformat_id, axis=1)

# Replace the 'position' column with the 'end' column
merged_df['position'] = merged_df['end']

# Remove the unnecessary columns
merged_df.drop(columns=['chr_x', 'start', 'end'], inplace=True)

# Reorder the columns to match the desired output
output_columns = ['id', 'chr', 'position', 'a0', 'a1', 'TYPE', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'ALL']
final_df = merged_df[output_columns]

# Save the reformatted output to a new file with space-separated values
final_df.to_csv('GRCh38_1000GP_Phase3_combined.legend', sep=' ', index=False)

print("Reformatting completed. The output file is 'GRCh38_1000GP_Phase3_combined.legend'.")

