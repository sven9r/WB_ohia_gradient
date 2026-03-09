# Gleason vs Clements Source file ---------------------------------------
# Setup script to load libraries, data, and define output paths
libs <- c("openxlsx","janitor", "vegan", "ggtext", 
          "RColorBrewer","glue", "here", "ggimage", "scales", 
          "viridis","ape", "phangorn", "parallel",
          "lme4", "patchwork", "ggsci", "ggvenn", "rgbif",
          "adespatial", "betapart", "geosphere", "ComplexUpset",
          "seqinr", "gt", "emmeans", "sf", "rnaturalearth" ,"rnaturalearthdata",
          "rnaturalearthhires",
          "ggspatial", "ggrepel", "ggnewscale", "elevatr", "raster","ggsci", "ggpubr",
          "classInt", "tidyverse", "ggh4x", "multcompView",
          "AICcmodavg", "segmented", "mgcv", "networkD3", "htmlwidgets")



# Install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}

# Load libraries
invisible(lapply(libs, library, character.only = T))


#WARNING
if (packageVersion("ggplot2") >= "4.0.0") {
  stop("ComplexUpset currently conflicts with ggplot2 4.0.0 during print(). Use ggplot2 3.5.2 for now.")
}

rm(installed_libs)
rm(libs)

# # Only resolve the conflicts you care about
# #conflicts
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("as.data.frame", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("union", "base")


#  Output directories and changes ----------------------------
# old outputs needed for legacy reasons

legacy_output <- here::here("output")   
legacy_figures <- here::here("figures") 
legacy_data  <- here::here("data")  

data_path <- function(...) file.path(legacy_data, ...)

# new output directories for the splitted papers --------------
out_root     <- here::here("papers", "gleason_clements_paper", "outputs")
out_fig_main <- file.path(out_root, "fig_main")
out_fig_supp <- file.path(out_root, "fig_supp")
out_tables   <- file.path(out_root, "tables")

dir.create(out_fig_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables,   recursive = TRUE, showWarnings = FALSE)


# Load the metadata ----------------
meta_data <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "metadata") |> 
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
  mutate(elevation_m = elevation_ft * 0.3048) |> 
  mutate(
    sample_site = dplyr::case_when(
      sample_site %in% c("Aiea", "'Aiea")            ~ "\u02bbAiea",   # ʻAiea
      sample_site %in% c("Ka'ala", "Kaala")          ~ paste0("Ka", "\u02bb", "ala"),  # Kaʻala
      sample_site %in% c("Schofield-Waikane", "Waikane")       ~ "Schoefield-Waik\u0101ne", # Waikāne
      TRUE ~ sample_site
    )
  ) 

meta_data_id <- meta_data |>
  dplyr::select(sample_id, sample_id_short) |> 
  distinct()

# Load zOTU table (ANML) ----------------
zotu_df <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "zotu_table_merged") |> 
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

# Define control pattern or list -----------
control_pattern <- "control|blank|mock|negative|water|buffer|reagent"

zotu_df_cleaned <- zotu_df_long |>
  mutate(sample_type = if_else(str_detect(sample_id_short, control_pattern), "control", "sample")) |>
  # Summarize max reads per zOTU in control samples
  group_by(zotu_id, sample_type) |> 
  reframe(control_reads = max(abundance)) |>
  # Keep only control reads
  filter(sample_type == "control") |> 
  dplyr::select(zotu_id, control_reads) |> 
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
  dplyr::select(sample_id_short, zotu_id, abundance_adjusted) |> 
  dplyr::filter(abundance_adjusted >= 5) |>
  pivot_wider(names_from = zotu_id, values_from = abundance_adjusted, values_fill = 0)



# Load the taxonomy table ----------
taxonomy_df <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "taxonomy") |> 
  janitor::clean_names() 

#alpha diversity ---------
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

# Join alpha diversity back to metadata ----------------
meta_data <- meta_data |> 
  left_join(alpha_div, by = "sample_id_short") |> 
  dplyr::select(-sample_id) |> 
  dplyr::rename(sample_id = sample_id_short) |> 
  dplyr::select(-sample_id_long) |> 
  distinct(sample_id, .keep_all = T)

meta_oahu <- meta_data %>%
  filter(island == "Oahu") |> 
  filter(!str_detect(sample_id, "CH|H")) |>
  filter(!str_detect(sample_id, "KT1")) |>
  filter(!str_detect(sample_id, "KA5")) |>
  filter(!str_detect(sample_id, "K_|G_"))

# Calculate jenks bands for elevation and rain --------------------
# Jenks natural breaks for 3 bands
jenks_brks <- classIntervals(meta_oahu$elevation_m, n = 3, style = "jenks")$brks
jenks_brks <- floor(jenks_brks) 
#add +1 to the last entry, buggy for the scale in plot
jenks_brks[length(jenks_brks)] <- jenks_brks[length(jenks_brks)] + 1

# --- Rainfall Jenks breaks (3 classes) ---
jenks_brks_rain <- classIntervals(meta_oahu$rainfall, n = 3, style = "jenks")$brks
jenks_brks_rain <- floor(jenks_brks_rain)

# add +1 to the last entry (avoid scale bug in plots)
jenks_brks_rain[length(jenks_brks_rain)] <- jenks_brks_rain[length(jenks_brks_rain)] + 1

# optional: wrap into a shorter alias for cutting
rain_brks <- jenks_brks_rain


meta_oahu <- meta_oahu %>%
  mutate(elev_band = cut(round(elevation_m, 0),
                         breaks = jenks_brks,
                         include.lowest = TRUE,
                         labels = c("Low", "Mid", "High"))) |> 
  mutate(rain_band = cut(round(rainfall, 0),
                         breaks = jenks_brks_rain,
                         include.lowest = TRUE,
                         labels = c("Low", "Mid", "High")))

# Optional: check band sizes
table(meta_oahu$elev_band)
table(meta_oahu$rain_band)


jenks_elevation_band <- meta_oahu |> 
  dplyr::select(sample_id,  elev_band) 

jenks_rain_band <- meta_oahu |> 
  dplyr::select(sample_id, rain_band) 

meta_data <- meta_data |>  
  left_join(jenks_elevation_band, by = "sample_id") 

zotu_df <- zotu_df_cleaned |> 
  dplyr::rename(sample_id = sample_id_short)


# helper for ggsave -----------------
fig_main_path <- function(name) file.path(out_fig_main, name)
fig_supp_path <- function(name) file.path(out_fig_supp, name)
tab_path      <- function(name) file.path(out_tables, name)

# Simple, canonical figure saver --------------------

save_fig <- function(
    tag,
    plot   = ggplot2::last_plot(),
    which  = c("main", "supp"),
    prefix = c("date", "none"),
    width  = 7,
    height = 5,
    bg = "white",
    dpi    = 300,
    ...
) {
  which  <- match.arg(which)
  prefix <- match.arg(prefix)
  
  dir <- if (which == "main") out_fig_main else out_fig_supp
  
  stamp <- if (prefix == "date") format(Sys.Date(), "%Y-%m-%d") else NULL
  filename <- paste(c(stamp, tag), collapse = "_")
  filename <- paste0(filename, ".png")
  
  ggplot2::ggsave(
    filename = file.path(dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg,
    ...
  )
}

# remove unused objects ------------------
rm(alpha_div,controls, control_pattern,zotu_df_long, meta_data_id,zotu_df_cleaned, jenks_brks, jenks_elevation_band, jenks_brks_rain, jenks_rain_band)
