import pandas as pd

# Define the chunk size for reading the file
chunksize = 5000  # Adjust based on your available memory

# Initialize an empty list to collect processed chunks
chunk_list = []

# Process the file in chunks
for chunk in pd.read_csv("../snp_array/TOPMED_imput/final_filtered_data.genotype.raw", sep=" ", header=None, chunksize=chunksize):

    # Set the IID as row names
    chunk.set_index(chunk[1], inplace=True)

    # Remove the first six columns (FID, IID, PAT, MAT, SEX, PHENOTYPE)
    chunk_data = chunk.iloc[:, 6:]

    # Append the processed chunk to the list
    chunk_list.append(chunk_data)

# Concatenate all chunks to form the full dataset
genotype_data = pd.concat(chunk_list, axis=0)

# Transpose the data to get SNPs as rows and samples as columns
genotype_data = genotype_data.T

# Set column names as IID
genotype_data.columns = genotype_data.index

# Save the formatted data to a file
genotype_data.to_csv("../snp_array/final_filtered_genotype_data.txt", sep=" ", header=True, index=False, quoting=False)

