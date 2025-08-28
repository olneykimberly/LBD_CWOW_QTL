import pandas as pd

# Read the merged output file
merged_df = pd.read_csv('merged_output.txt', sep='\t')

# Rename columns
merged_df.rename(columns={'chr_x': 'chr'}, inplace=True)

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
merged_df.drop(columns=['chr_y', 'start', 'end'], inplace=True)

# Reorder the columns to match the desired output
output_columns = ['id', 'chr', 'position', 'a0', 'a1', 'TYPE', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'ALL']
final_df = merged_df[output_columns]

# Save the reformatted output to a new file with space-separated values
final_df.to_csv('reformatted_output.txt', sep=' ', index=False)

print("Reformatting completed. The output file is 'reformatted_output.txt'.")

