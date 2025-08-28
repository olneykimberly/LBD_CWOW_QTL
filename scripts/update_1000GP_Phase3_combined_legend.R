# Load necessary libraries
library(data.table)
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/QTL/LBD_CWOW_QTL/snp_array/reference/eagle_1k_reference")
# Read the files
legend <- fread("1000GP_Phase3_combined.legend", header = TRUE)
lifted <- fread("output.bed", header = FALSE, col.names = c("chr", "new_position_start", "new_position_end", "id"))

# Split the ID into components
lifted[, `:=`(original_position = sub(".*:", "", sub(":.*", "", id)),
              a0 = sub(".*:(.*?):.*", "\\1", id),
              a1 = sub(".*:.*?:(.*)", "\\1", id))]

# Merge the legend data with the lifted positions based on the id
merged_data <- merge(lifted, legend, by.x = "id", by.y = "id")

# Update the ID with the new position
merged_data[, new_id := paste0(chr.x, ":", new_position_end, ":", a1.x)]

# Reorder and format the final output
final_output <- merged_data[, .(id = new_id, chr = chr.x, position = new_position_end, a0 = a0.x, a1 = a1.x, TYPE, AFR, AMR, EAS, EUR, SAS, ALL)]


# Remove values after the second : and before the third : in the id column
final_output[, id := sub("([^:]*:[^:]*):[^:]*:(.*)", "\\1:\\2", id)]
head(final_output)
# Remove the values in the a0 column
final_output[, a0 := NULL]

# Assuming final_output is already a data.table
setDT(final_output)

# Split the a1 column into two new columns a0 and a1
# Split the a1 column into a list of vectors
split_values <- tstrsplit(final_output$a1, ":", fixed = TRUE)
# Convert the list into a data.table with appropriate column names
split_dt <- data.table(do.call(rbind, split_values))

# Rename columns based on the number of splits
setnames(split_dt, c("a0", "a1"))

head(final_output)

# Write the final output to a file
fwrite(final_output, "updated_legend_GRCh38.txt", sep = " ", col.names = TRUE)
