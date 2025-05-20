###########################################################
#
#   Title:  
#   ohia gradient
# 
#   Source file for the repository
#   soon to be found here
#
#   Author: 
#   Sven Weber
#
#   Start date:   
#   02/2025
#
#   version info:
#   base R v: 4.4.2 | 
#   Rstudio v: 2025.2.1 |
#   Rmarkdown v: 2.11 |
#   Renv v: NOT RUNNING
#   
############################################################



# Load necessary libraries
libs <- c("tidyverse","openxlsx","janitor", 
          "vegan", "ggtext", "RColorBrewer",
          "glue", "here", "ggimage", "scales", 
          "viridis","Biostrings", "ape", 
          "phangorn", "parallel")

# Load libraries
invisible(lapply(libs, library, character.only = T))

rm(libs)

# Load the metadata
metadata <- openxlsx::read.xlsx(here::here("data/raw_data/ANML.xlsx"), sheet = "meta_data" ) |> 
  dplyr::select(1:13) 
  

# Load the taxonomy table
taxonomy_unfiltered <- openxlsx::read.xlsx(here::here("data/raw_data/ANML.xlsx"), sheet = "tax_unfiltered") 

# Count occurrences across all columns
taxonomy_counts <- taxonomy_unfiltered |> 
  group_by(query, perc_id, kingdom, phylum, superclass, class, subclass, 
          order, infraorder, superfamily, family, genus, species) %>%
  reframe(n = n()) |> 
  select(query, perc_id, n, everything())

tax_filter = 91

best_taxonomy <- taxonomy_counts |> 
  filter(perc_id >= tax_filter) |>
  group_by(query) |> 
  arrange(desc(perc_id)) |>
  mutate(id_diff = perc_id - lead(perc_id, default = 0),.before =3,
         same_family = ifelse(family == lead(family), "yes", lead(family, default = "NA")),
         same_genus = ifelse(genus == lead(genus), "yes", lead(genus, default = "NA")),
         same_species = ifelse(species == lead(species), "yes", lead(species, default = "NA"))) |>
  slice_max(n = 1, order_by = across(c(perc_id, n)), with_ties = F) |> 
  ungroup() |> 
  filter(phylum == "Arthropoda") |> 
  mutate(zotu_num = as.integer(str_replace(query, "Zotu", "")),.before =2) |>
  arrange(zotu_num) |>
  dplyr::rename(zotu_id = "query")

best_tax_short <- best_taxonomy |>
  select(zotu_id)

# Load zOTU table (ANML)
# Merge with filtered taxonomy table to reduce to arthropods reads
zotu_df_unfiltered <- openxlsx::read.xlsx(here::here("data/raw_data/ANML.xlsx"), sheet = "otu_table_unfiltered") |>
  pivot_longer(cols = -sample_id_long, names_to = "zotu_id", values_to = "reads") |>
  inner_join(best_tax_short, by = "zotu_id") |>
  pivot_wider(names_from = zotu_id, values_from = reads) 


zotu_df_merged <- zotu_df_unfiltered |> 
  pivot_longer(cols = -sample_id_long, names_to = "zotu_id", values_to = "reads") |>
  left_join(metadata, by = "sample_id_long") |>
  group_by(sample_id, zotu_id) |> 
  reframe(reads = sum(reads)) |>
  pivot_wider(names_from = zotu_id, values_from = reads) |>
  select(sample_id, starts_with("Zotu")) 



#write new files in the output folder so we can access and preload everything according to our needs from one single source file
#also make sure the output is dynamic and reflects at least the day when you ran this preprocessing file
#otherwise you might end up with outdated files and or overwrite files from prior runs


#saves current date into a variable
current_date <- format(Sys.Date(), "%Y-%m-%d")
#use glue to create a dynamic filename
filename_tax <- glue("output/{current_date}_tax_table_filtered_{current_date}.csv")
filename_df_single <- glue("output/{current_date}_otu_table_filtered_{current_date}.csv")
filename_df_merged <- glue("output/{current_date}_otu_table_filtered_{current_date}.csv")
filename_meta_table <- glue("output/{current_date}_meta_table_{current_date}.csv")

# combine created files into .xlsx for convince but keep dated files for reproducibility
#save the new and filtered data frames in one file with openxlsx
wb <- createWorkbook()
addWorksheet(wb, "zotu_table_single")
addWorksheet(wb, "zotu_table_merged")
addWorksheet(wb, "metadata")
addWorksheet(wb, "taxonomy")
writeData(wb, "taxonomy", best_taxonomy)
writeData(wb, "zotu_table_single", zotu_df_unfiltered)
writeData(wb, "zotu_table_merged", zotu_df_merged)
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file = glue("data/ohia_gradient_{current_date}_{tax_filter}.xlsx"), overwrite = F)
