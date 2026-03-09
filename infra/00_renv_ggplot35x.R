source("scripts/00_ggplotforce.R")


# 00_install_all_after_pin.R
# Pin ggplot2 = 3.5.2, ggimage = 0.3.3, then install ALL requested libs.

# --- your full library list (alphabetical, deduped) ---
libs <- sort(unique(c(
  "adespatial","ape","betapart","classInt","ComplexUpset",
  "elevatr","emmeans","geosphere","geepack","ggimage","ggnewscale",
  "ggpubr","ggrepel","ggh4x","ggsci","ggspatial","ggtext",
  "ggvenn","glue","gt","here","janitor","lme4","openxlsx",
  "parallel","patchwork","phangorn","raster","RColorBrewer",
  "rgbif","rnaturalearth","rnaturalearthdata","scales",
  "seqinr","sf","tidyverse","vegan","viridis"
)))

# Helper: install a specific version from CRAN Archive
.install_arch <- function(pkg, ver) {
  url <- sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", pkg, pkg, ver)
  message(">> Installing ", pkg, " ", ver, " from: ", url)
  renv::install(url)
}

# 0) Activate renv (init if first time)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv.lock")) renv::init(bare = TRUE) else renv::activate()

# 1) Pin ggplot2 3.5.2 (+ core companions) and restore clean
renv::record(list(
  ggplot2 = "3.5.2",
  scales  = "1.3.0",
  gtable  = "0.3.5",
  isoband = "0.2.7",
  rlang   = "1.1.4",
  vctrs   = "0.6.5",
  cli     = "3.6.3",
  glue    = "1.7.0",
  tibble  = "3.2.1",
  withr   = "3.0.0"
))
renv::restore(prompt = FALSE, clean = TRUE)

# 2) ggimage: force pre-4.0.0-compatible version
.install_arch("ggimage", "0.3.3")
renv::record(list(ggimage = as.character(packageVersion("ggimage"))))

# 3) Install the rest of your packages (excluding ggimage we just installed)
remaining <- setdiff(libs, "ggimage")
renv::install(remaining)

# 4) Snapshot final environment
renv::snapshot(prompt = FALSE)

# 5) Quick sanity print
cat("\nPinned:\n",
    "  ggplot2  ", as.character(packageVersion("ggplot2")), "\n",
    "  ggimage  ", as.character(packageVersion("ggimage")),  "\n\n",
    "Sample checks:\n",
    "  ggh4x    ", if (requireNamespace("ggh4x", quietly = TRUE)) as.character(packageVersion("ggh4x")) else "not installed", "\n",
    "  ggpubr   ", if (requireNamespace("ggpubr", quietly = TRUE)) as.character(packageVersion("ggpubr")) else "not installed", "\n",
    "  ComplexUpset ", if (requireNamespace("ComplexUpset", quietly = TRUE)) as.character(packageVersion("ComplexUpset")) else "not installed", "\n",
    sep = "")

