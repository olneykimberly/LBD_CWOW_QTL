import pandas as pd

# Read the input file
input_file = 'updated_legend_GRCh38.txt'  # Replace with your input file path
df = pd.read_csv(input_file, sep=' ')  # Adjust 'sep' if needed

# Strip whitespace from column names
df.columns = df.columns.str.strip()

# Print column names to verify
print("Column names:", df.columns)

# Check if 'id' and 'a1' columns are present
if 'id' not in df.columns or 'a1' not in df.columns:
    raise KeyError("'id' or 'a1' column is missing from the DataFrame")

# Remove the value between the second and third ":" in the 'id' column
df['id'] = df['id'].apply(lambda x: ':'.join(x.split(':')[:2] + x.split(':')[3:]))

# Drop the 'a0' column
if 'a0' in df.columns:
    df = df.drop(columns=['a0'])

# Ensure 'a1' column is of string type
df['a1'] = df['a1'].astype(str)

# Split the 'a1' column by ":" and create temporary columns
split_columns = df['a1'].str.split(':', expand=True)

# Filter out rows where the split does not return exactly 2 columns
valid_split = split_columns.shape[1] == 2
if valid_split:
    # Keep only rows where split resulted in exactly 2 columns
    valid_rows = split_columns[split_columns.notna().all(axis=1)]
    
    # Assign the split columns back to the DataFrame
    df = df.loc[valid_rows.index]
    df[['a0', 'a1']] = valid_rows

# Drop the original 'a1' column
df = df.drop(columns=['a1'])

# Save the updated DataFrame to a new file
output_file = 'output_file.txt'  # Replace with your desired output file path
df.to_csv(output_file, sep=' ', index=False)  # Adjust 'sep' if needed

