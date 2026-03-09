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
          "viridis", "ape", 
          "phangorn", "parallel")

# Load libraries
invisible(lapply(libs, library, character.only = T))

rm(libs)

# Load the metadata
metadata <- openxlsx::read.xlsx(here::here("data/raw_data/ANML.xlsx"), sheet = "meta_data" ) |> 
  dplyr::select(1:14) 
  

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


taxonomy_BOLD <- openxlsx::read.xlsx(here::here("data/raw_data/boldigger3_data/ANML_BOLD_result.xlsx")) |> 
  janitor::clean_names() |>
  filter(phylum == "Arthropoda") |> 
  mutate(otu_num = as.numeric(str_replace(id, "Zotu", "")),.before =2) |>
  arrange(otu_num) |> 
  as.data.frame()

# BOLD VS NCBI comparison
# filter best_taxonomy the same way BOLDIGGER standard works
# Species: 99, Genus: 97, Family: 90, Order: 85 fill information above with blanks when perc_id is below the threshold
# create a new column that indicates which level of taxonomy we can trust based on the perc_id for the NCBI data

# --- hard prefilter for comparison logic ---
tax_prefilter <- 84

prep_comparison <- taxonomy_counts |>
  filter(perc_id >= tax_prefilter) |>
  group_by(query) |>
  arrange(desc(perc_id), desc(n)) |>
  #remove all columns with "species abbreviations" to avoid issues with the following logics with str_detect
   mutate(across(c(kingdom, phylum, superclass, class, subclass, order, infraorder, superfamily, family, genus, species), ~if_else(str_detect(., " sp\\."), NA_character_, .))) |>
  #remove "Nan" from all columns to avoid issues with the following logics
   mutate(across(c(kingdom, phylum, superclass, class, subclass, order, infraorder, superfamily, family, genus, species), ~na_if(., "Nan"))) |>
  
  # decide selected_level based on best available hit PER query AND IF COLUMN ISNT EMPTY
  mutate(
    selected_level = case_when(
      perc_id >= 99 & !is.na(species) ~ "species",
      perc_id >= 97 & !is.na(genus)   ~ "genus",
      perc_id >= 90 & !is.na(family)  ~ "family",
      perc_id >= 85 & !is.na(order)   ~ "order",
      TRUE ~ "higher"
    )
  ) |>
  
  # ambiguity metrics computed BEFORE slicing
  mutate(
    id_diff = perc_id - lead(perc_id),
    same_family  = if_else(family  == lead(family),  "yes", "no", missing = "no"),
    same_genus   = if_else(genus   == lead(genus),   "yes", "no", missing = "no"),
    same_species = if_else(species == lead(species), "yes", "no", missing = "no"),
    .before = 3
  ) |>
  
  # now safely keep the best hit
  slice(1) |>
  ungroup() |>
  filter(phylum == "Arthropoda") |>
  mutate(zotu_num = as.integer(str_replace(query, "Zotu", "")), .before = 2) |>
  arrange(zotu_num) |>
  rename(zotu_id = query)


# final taxonomy table with epistemic masking
rank_order <- c(
  "kingdom", "phylum", "superclass", "class", "subclass",
  "order", "infraorder", "superfamily", "family", "genus", "species"
)

taxonomy_comparison <- prep_comparison |>
  mutate(
    across(
      all_of(rank_order),
      ~ if_else(
        match(cur_column(), rank_order) <= match(selected_level, rank_order),
        as.character(.x),
        NA_character_
      )
    )
  ) |>
  select(
    zotu_id, perc_id,
    all_of(rank_order),
    selected_level,
    id_diff, same_family, same_genus, same_species
  )


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
filename_tax_BOLD <- glue("output/{current_date}_tax_table_filtered2_{current_date}.xlsx")
filename_df_single <- glue("output/{current_date}_otu_table_filtered_{current_date}.csv")
filename_df_merged <- glue("output/{current_date}_otu_table_filtered_{current_date}.csv")
filename_meta_table <- glue("output/{current_date}_meta_table_{current_date}.csv")
filename_comparison <- glue("output/{current_date}_tax_comparison_{current_date}.csv")

# combine created files into .xlsx for convince but keep dated files for reproducibility
#save the new and filtered data frames in one file with openxlsx
wb <- createWorkbook()
addWorksheet(wb, "zotu_table_single")
addWorksheet(wb, "zotu_table_merged")
addWorksheet(wb, "metadata")
addWorksheet(wb, "taxonomy")
addWorksheet(wb, "taxonomy_BOLD")
addWorksheet(wb, "taxonomy_comparison")
writeData(wb, "taxonomy_comparison", taxonomy_comparison)
writeData(wb, "taxonomy", best_taxonomy)
writeData(wb, "taxonomy_BOLD", taxonomy_BOLD)
writeData(wb, "zotu_table_single", zotu_df_unfiltered)
writeData(wb, "zotu_table_merged", zotu_df_merged)
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file = glue("data/ohia_gradient_{current_date}.xlsx"), overwrite = F)


