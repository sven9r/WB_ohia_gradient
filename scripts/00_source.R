#libraries
libs <- c("openxlsx","janitor", "vegan", "ggtext", 
          "RColorBrewer","glue", "here", "ggimage", "scales", 
          "viridis","Biostrings", "ape", "phangorn", "parallel",
           "lme4", "patchwork", "ggsci", "ggvenn", "rgbif",
          "adespatial", "betapart", "geosphere", "ComplexUpset",
          "seqinr", "gt", "emmeans", "sf", "rnaturalearth" ,"rnaturalearthdata",
          "ggspatial", "ggrepel", "ggnewscale", "elevatr", "raster","ggsci", "ggpubr",
          "tidyverse")


# Install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}

# Load libraries
invisible(lapply(libs, library, character.only = T))

rm(installed_libs)
rm(libs)

# # Only resolve the conflicts you care about
# #conflicts
# library(conflicted)
# conflict_prefer("filter", "dplyr")
# conflict_prefer("select", "dplyr")
# conflict_prefer("as.data.frame", "base")


# Load the metadata
meta_data <- openxlsx::read.xlsx(here::here("data/ohia_gradient_2025-03-25_91.xlsx"), sheet = "metadata") |> 
  janitor::clean_names() |> 
  mutate(sample_id_short = sample_id %>%
           # Step 1: Replace tree codes with a single letter
           str_replace("^GA_", "G_") %>%    # GA -> G
           str_replace("^KO_", "K_") %>%    # Ko -> K
           str_replace("^OH_", "O_") %>%    # OH -> O
           # 2. Remove trailing digits in the last portion (e.g., AL1 → AL)
           str_replace("(_[A-Za-z]+)\\d+$", "\\1"),
         # 3. Add elevation code (1, 2, or 3) based on the middle portion of the sample ID
         elev_num = case_when(
           str_detect(sample_id_short,"_L_") ~ "1",
           str_detect(sample_id_short,"_M_") ~ "2",
           str_detect(sample_id_short,"_H_") ~ "3",
           TRUE ~ NA_character_
         )) |> 
  # 4. Remove the middle elevation code (_L_, _M_, or _H_)
  mutate(sample_id_short = sample_id_short |> 
           str_replace("_[LMH]_", "_")) |> 
  mutate(sample_id_short_new = if_else(
    is.na(elev_num),
    sample_id_short,
    paste0(sample_id_short, elev_num)
  ),.before = sample_id_short) |> 
  dplyr::select(-sample_id_short, -elev_num) |> 
  dplyr::rename(sample_id_short = sample_id_short_new) |>
  dplyr::select(sample_id, sample_id_short, everything()) |> 
  #calculate elevation_m instead of ft for publication
  mutate(elevation_m = elevation_ft * 0.3048) 

meta_data_id <- meta_data |>
  select(sample_id, sample_id_short) |> 
  distinct()

# Load zOTU table (ANML)
zotu_df <- openxlsx::read.xlsx(here::here("data/ohia_gradient_2025-03-25_91.xlsx"), sheet = "zotu_table_merged") |> 
  janitor::clean_names() 

zotu_df_long <- zotu_df |>
  pivot_longer(-sample_id, names_to = "zotu_id", values_to = "abundance") |> 
  #join the new sample_id 
  left_join(meta_data_id, by = "sample_id") |>
  group_by(sample_id_short, zotu_id) |>
  #calculate summed reads per replicates in Leke's data
  reframe(abundance = sum(abundance))

controls <- zotu_df |> 
  filter(str_detect(sample_id, "control"))

# Define control pattern or list
control_pattern <- "control|blank|mock|negative|water|buffer|reagent"

zotu_df_cleaned <- zotu_df_long |>
  mutate(sample_type = if_else(str_detect(sample_id_short, control_pattern), "control", "sample")) |>
  # Summarize max reads per zOTU in control samples
  group_by(zotu_id, sample_type) |> 
  reframe(control_reads = max(abundance)) |>
  # Keep only control reads
  filter(sample_type == "control") |> 
  select(zotu_id, control_reads) |> 
  # Join back to original data
  right_join(zotu_df_long, by = "zotu_id") |> 
  #get rid of zotu1, probably whole specimen
  filter(zotu_id != "zotu1") |>
  # Replace NA control reads with 0 (for OTUs not in any control)
  mutate(control_reads = replace_na(control_reads, 0)) |> 
  # Subtract control reads only for real samples
  mutate(abundance_adjusted = case_when(
    !str_detect(sample_id_short, control_pattern) ~ pmax(abundance - control_reads, 0),
    TRUE ~ abundance  # leave control reads unchanged
  )) |> 
  select(sample_id_short, zotu_id, abundance_adjusted) |> 
  filter(abundance_adjusted >= 5) |>
  pivot_wider(names_from = zotu_id, values_from = abundance_adjusted, values_fill = 0)



# Load the taxonomy table
taxonomy_df <- openxlsx::read.xlsx(here::here("data/ohia_gradient_2025-03-25_91.xlsx"), sheet = "taxonomy") |> 
  janitor::clean_names() 

#alpha diversity
alpha_div <- zotu_df_cleaned |> 
  column_to_rownames("sample_id_short") |>
  as.matrix() |> 
  (\(x) data.frame(
    shannon = vegan::diversity(x,index = "shannon"),  # Shannon diversity
    simpson = vegan::diversity(x,index = "simpson"),  # Simpson diversity
    richness = vegan::specnumber(x)  # Species richness (Observed ASVs)
  ))() |> 
  rownames_to_column(var = "sample_id_short") |>  # Add sample names back
  as.data.frame()


meta_data <- meta_data |> 
  left_join(alpha_div, by = "sample_id_short") |> 
  select(-sample_id) |> 
  dplyr::rename(sample_id = sample_id_short) |> 
  dplyr::select(-sample_id_long) |> 
  distinct(sample_id, .keep_all = T)
  
zotu_df <- zotu_df_cleaned |> 
  dplyr::rename(sample_id = sample_id_short)


rm(alpha_div,controls, control_pattern,zotu_df_long, meta_data_id,zotu_df_cleaned)
