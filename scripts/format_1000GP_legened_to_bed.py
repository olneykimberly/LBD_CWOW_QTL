import pandas as pd

# Load the files
legend = pd.read_csv("1000GP_Phase3_combined.legend", delim_whitespace=True)
lifted = pd.read_csv("output.bed", sep='\t', header=None, names=["chr", "new_position_start", "new_position_end", "id"])

# Extract information from the 'id' column in lifted
lifted[['original_chr', 'original_position', 'a0', 'a1']] = lifted['id'].str.split(':', expand=True)

# Merge based on the 'id' from legend and lifted
merged = pd.merge(lifted, legend, left_on='id', right_on='id', how='inner')

# Create new 'id' with updated position
merged['new_id'] = merged['chr'].astype(str) + ":" + merged['new_position_end'].astype(str) + ":" + merged['a0'] + ":" + merged['a1']

# Reorder the columns to match the legend format
final_output = merged[['new_id', 'chr', 'new_position_end', 'a0', 'a1', 'TYPE', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'ALL']]

# Rename the columns to match the legend format
final_output.columns = ['id', 'chr', 'position', 'a0', 'a1', 'TYPE', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'ALL']

# Save the output to a new file
final_output.to_csv("updated_legend_GRCh38.txt", sep=' ', index=False)

