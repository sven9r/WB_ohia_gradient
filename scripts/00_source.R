# Project analysis source ---------------------------------------

libs <- c(
  "openxlsx", "janitor", "vegan", "ggtext", "RColorBrewer", "glue",
  "here", "ggimage", "scales", "viridis",
  "lme4", "patchwork", "ggsci", "ggvenn",
  "betapart", "geosphere", "ComplexUpset", "gt", "emmeans",
  "sf", "ggspatial", "ggrepel", "ggnewscale", "raster",
  "ggpubr", "classInt", "tidyverse", "ggh4x", "multcompView",
  "networkD3", "htmlwidgets"
)

missing_libs <- libs[!vapply(libs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_libs)) {
  stop("Missing packages: ", paste(missing_libs, collapse = ", "), ". Run renv::restore() from the repository root.")
}

invisible(lapply(libs, library, character.only = TRUE))

if (packageVersion("ggplot2") >= "4.0.0") {
  stop("ComplexUpset currently conflicts with ggplot2 4.0.0 during print(). Use ggplot2 3.5.2 for now.")
}

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("as.data.frame", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("union", "base")

rm(missing_libs)
rm(libs)

legacy_output <- here::here("output")
legacy_figures <- here::here("figures")
legacy_data <- here::here("data")

data_path <- function(...) file.path(legacy_data, ...)

out_root <- legacy_figures
out_fig_main <- file.path(out_root, "main")
out_fig_supp <- file.path(out_root, "supp")
out_fig_exploratory <- file.path(out_root, "exploratory")
out_tables <- legacy_output

dir.create(out_fig_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_exploratory, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

meta_data <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "metadata") |>
  janitor::clean_names() |>
  mutate(
    sample_id_short = sample_id %>%
      str_replace("^GA_", "G_") %>%
      str_replace("^KO_", "K_") %>%
      str_replace("^OH_", "O_") %>%
      str_replace("(_[A-Za-z]+)\\d+$", "\\1"),
    elev_num = case_when(
      str_detect(sample_id_short, "_L_") ~ "1",
      str_detect(sample_id_short, "_M_") ~ "2",
      str_detect(sample_id_short, "_H_") ~ "3",
      TRUE ~ NA_character_
    )
  ) |>
  mutate(sample_id_short = sample_id_short |> str_replace("_[LMH]_", "_")) |>
  mutate(
    sample_id_short_new = if_else(
      is.na(elev_num),
      sample_id_short,
      paste0(sample_id_short, elev_num)
    ),
    .before = sample_id_short
  ) |>
  dplyr::select(-sample_id_short, -elev_num) |>
  dplyr::rename(sample_id_short = sample_id_short_new) |>
  dplyr::select(sample_id, sample_id_short, everything()) |>
  mutate(
    elevation_m = elevation_ft * 0.3048,
    sample_site = dplyr::case_when(
      sample_site %in% c("Aiea", "'Aiea") ~ "\u02bbAiea",
      sample_site %in% c("Ka'ala", "Kaala") ~ paste0("Ka", "\u02bb", "ala"),
      sample_site %in% c("Manana", "Mānana") ~ "Mānana",
      sample_site %in% c("Schofield-Waikane", "Waikane") ~ "Schofield-Waik\u0101ne",
      TRUE ~ sample_site
    )
  )

meta_data_id <- meta_data |>
  dplyr::select(sample_id, sample_id_short) |>
  distinct()

zotu_df <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "zotu_table_merged") |>
  janitor::clean_names()

zotu_df_long <- zotu_df |>
  pivot_longer(-sample_id, names_to = "zotu_id", values_to = "abundance") |>
  left_join(meta_data_id, by = "sample_id") |>
  group_by(sample_id_short, zotu_id) |>
  reframe(abundance = sum(abundance))

controls <- zotu_df |>
  filter(str_detect(sample_id, "control"))

control_pattern <- "control|blank|mock|negative|water|buffer|reagent"

zotu_df_cleaned <- zotu_df_long |>
  mutate(sample_type = if_else(str_detect(sample_id_short, control_pattern), "control", "sample")) |>
  group_by(zotu_id, sample_type) |>
  reframe(control_reads = max(abundance)) |>
  filter(sample_type == "control") |>
  dplyr::select(zotu_id, control_reads) |>
  right_join(zotu_df_long, by = "zotu_id") |>
  filter(zotu_id != "zotu1") |>
  mutate(control_reads = replace_na(control_reads, 0)) |>
  mutate(
    abundance_adjusted = case_when(
      !str_detect(sample_id_short, control_pattern) ~ pmax(abundance - control_reads, 0),
      TRUE ~ abundance
    )
  ) |>
  dplyr::select(sample_id_short, zotu_id, abundance_adjusted) |>
  dplyr::filter(abundance_adjusted >= 5) |>
  pivot_wider(names_from = zotu_id, values_from = abundance_adjusted, values_fill = 0)

taxonomy_df <- openxlsx::read.xlsx(data_path("ohia_gradient_2026-02-09.xlsx"), sheet = "taxonomy") |>
  janitor::clean_names()

alpha_div <- zotu_df_cleaned |>
  column_to_rownames("sample_id_short") |>
  as.matrix() |>
  (\(x) data.frame(
    shannon = vegan::diversity(x, index = "shannon"),
    simpson = vegan::diversity(x, index = "simpson"),
    richness = vegan::specnumber(x)
  ))() |>
  rownames_to_column(var = "sample_id_short") |>
  as.data.frame()

meta_data <- meta_data |>
  left_join(alpha_div, by = "sample_id_short") |>
  dplyr::select(-sample_id) |>
  dplyr::rename(sample_id = sample_id_short) |>
  dplyr::select(-sample_id_long) |>
  distinct(sample_id, .keep_all = TRUE)

meta_oahu <- meta_data %>%
  filter(island == "Oahu") |>
  filter(!str_detect(sample_id, "CH|H")) |>
  filter(!str_detect(sample_id, "KT1")) |>
  filter(!str_detect(sample_id, "KA5")) |>
  filter(!str_detect(sample_id, "K_|G_"))

jenks_brks <- classIntervals(meta_oahu$elevation_m, n = 3, style = "jenks")$brks
jenks_brks <- floor(jenks_brks)
jenks_brks[length(jenks_brks)] <- jenks_brks[length(jenks_brks)] + 1

jenks_brks_rain <- classIntervals(meta_oahu$rainfall, n = 3, style = "jenks")$brks
jenks_brks_rain <- floor(jenks_brks_rain)
jenks_brks_rain[length(jenks_brks_rain)] <- jenks_brks_rain[length(jenks_brks_rain)] + 1

rain_brks <- jenks_brks_rain

meta_oahu <- meta_oahu %>%
  mutate(
    elev_band = cut(
      round(elevation_m, 0),
      breaks = jenks_brks,
      include.lowest = TRUE,
      labels = c("Low", "Mid", "High")
    )
  ) |>
  mutate(
    rain_band = cut(
      round(rainfall, 0),
      breaks = jenks_brks_rain,
      include.lowest = TRUE,
      labels = c("Low", "Mid", "High")
    )
  )

jenks_elevation_band <- meta_oahu |>
  dplyr::select(sample_id, elev_band)

jenks_rain_band <- meta_oahu |>
  dplyr::select(sample_id, rain_band)

meta_data <- meta_data |>
  left_join(jenks_elevation_band, by = "sample_id")

zotu_df <- zotu_df_cleaned |>
  dplyr::rename(sample_id = sample_id_short)

fig_main_path <- function(name) file.path(out_fig_main, name)
fig_supp_path <- function(name) file.path(out_fig_supp, name)
fig_exploratory_path <- function(name) file.path(out_fig_exploratory, name)
tab_path <- function(name) file.path(out_tables, name)

save_fig <- function(
    tag,
    plot = ggplot2::last_plot(),
    which = c("main", "supp", "exploratory"),
    prefix = c("date", "none"),
    width = 7,
    height = 5,
    bg = "white",
    dpi = 300,
    ...
) {
  which <- match.arg(which)
  prefix <- match.arg(prefix)

  date_token <- format(Sys.Date(), "%Y_%m_%d")
  date_parts <- strsplit(date_token, "_", fixed = TRUE)[[1]]

  parts <- strsplit(tag, "_", fixed = TRUE)[[1]]

  # Support explicit placeholder tags such as "F1_date_map".
  if ("date" %in% parts) {
    expanded <- character(0)
    for (p in parts) {
      if (identical(p, "date")) {
        expanded <- c(expanded, date_parts)
      } else {
        expanded <- c(expanded, p)
      }
    }
    parts <- expanded
  }

  has_ymd_triplet <- FALSE
  if (length(parts) >= 3) {
    for (i in seq_len(length(parts) - 2)) {
      if (
        grepl("^[0-9]{4}$", parts[i]) &&
        grepl("^[0-9]{2}$", parts[i + 1]) &&
        grepl("^[0-9]{2}$", parts[i + 2])
      ) {
        has_ymd_triplet <- TRUE
        break
      }
    }
  }

  # Ensure recursive date insertion for all figure tags.
  if (!has_ymd_triplet) {
    parts <- append(parts, date_parts, after = 1)
  }

  tag_out <- paste(parts, collapse = "_")

  dir <- dplyr::case_match(
    which,
    "main" ~ out_fig_main,
    "supp" ~ out_fig_supp,
    "exploratory" ~ out_fig_exploratory
  )

  stamp <- if (prefix == "date" && !grepl("[0-9]{4}_[0-9]{2}_[0-9]{2}", tag_out)) date_token else NULL
  filename <- paste(c(stamp, tag_out), collapse = "_")
  filename <- paste0(filename, ".png")

  ggplot2::ggsave(
    filename = file.path(dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    ...
  )
}

rm(
  alpha_div, controls, control_pattern, zotu_df_long, meta_data_id,
  zotu_df_cleaned, jenks_brks, jenks_elevation_band, jenks_brks_rain,
  jenks_rain_band
)
