# Required libraries
library(data.table)

# Define input and output directories
input_dir <- "/research/labs/neurology/fryer/m239830/LBD_CWOW/QTL/LBD_CWOW_QTL/MatrixEQTL"  
output_dir <- "/research/labs/neurology/fryer/m239830/LBD_CWOW/QTL/shiny_app"  

# List of filenames to process
file_names <- c("cis_eQTL_LBD_control", "trans_eQTL_LBD_control",
                "cis_eQTL_LBD_ATS_control", "trans_eQTL_LBD_ATS_control",
                "cis_eQTL_LBD_AS_control", "trans_eQTL_LBD_AS_control",
                "cis_eQTL_LBD_S_control", "trans_eQTL_LBD_S_control",
                "cis_eQTL_AD_control", "trans_eQTL_AD_control",
                "cis_eQTL_PA_control", "trans_eQTL_PA_control")

# FDR threshold
fdr_threshold <- 0.001

# Loop through each file, filter, and save
for (file_name in file_names) {
  
  # Construct full file path
  input_file_path <- file.path(input_dir, paste0(file_name))
  
  # Check if file exists
  if (file.exists(input_file_path)) {
    
    # Read the data
    data <- fread(input_file_path)
    
    # Filter rows by FDR threshold
    filtered_data <- data[FDR < fdr_threshold]
    
    # Construct output file path
    output_file_path <- file.path(output_dir, paste0(file_name))
    
    # Save the filtered data
    write.table(filtered_data, output_file_path, quote = FALSE, row.names = FALSE, sep ="\t")
    
    cat("Filtered data saved to:", output_file_path, "\n")
    
  } else {
    cat("File not found:", input_file_path, "\n")
  }
}
