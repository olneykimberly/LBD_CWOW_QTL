map <- fread("../snp_array/reference/CHRPOSREFALT_to_rsid.txt", header = FALSE)
setnames(map, c("oldID","rsid"))

# remove rows without rsID
map <- map[!is.na(rsid)]

fwrite(
  map,
  "../snp_array/reference/update_rsid_map.txt",
  sep = "\t",
  col.names = FALSE
)
