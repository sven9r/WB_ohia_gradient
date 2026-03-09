# who_forces_ggplot4.R
# Find which package(s) force ggplot2 >= 4.0.0 (directly or via deps).
# Uses CRAN metadata only; does NOT install anything.

suspects <- c(
  "ggimage","ggnewscale","ggpubr","ggrepel","ggh4x",
  "ggsci","ggspatial","ggtext","ggvenn"
)

target <- "3.5.1"

# Use your current CRAN; if you prefer, freeze a snapshot here:
# options(repos = c(CRAN = "https://cran.r-project.org"))

ap <- utils::available.packages()  # CRAN index

# Pull the dependency strings for a package from available.packages()
dep_strings <- function(pkg) {
  if (!pkg %in% rownames(ap)) return(character())
  fields <- ap[pkg, c("Depends","Imports","LinkingTo","Suggests")]
  fields[is.na(fields)] <- ""
  fields
}

# Extract all (pkg, constraint) pairs from a DESCRIPTION-style dep string
parse_deps <- function(x) {
  x <- paste(x, collapse = ",")
  x <- gsub("\\n", " ", x)
  parts <- unlist(strsplit(x, ","))
  parts <- trimws(parts[nzchar(parts)])
  # split "name (>= 1.2.3)" into name + constraint
  name  <- sub("\\s*\\(.*\\)$", "", parts)
  const <- ifelse(grepl("\\(", parts),
                  sub(".*\\((.*)\\).*", "\\1", parts),
                  "")
  data.frame(pkg = name, constraint = const, stringsAsFactors = FALSE)
}

# Recursively gather all hard deps for a package
all_deps <- function(pkg) {
  # recursive graph from tools (Depends/Imports/LinkingTo; no Suggests)
  g <- tools::package_dependencies(packages = pkg, db = ap, which = c("Depends","Imports","LinkingTo"), recursive = TRUE)
  unique(unlist(g, use.names = FALSE))
}

# Pull the strongest ggplot2 requirement seen across a package + its deps
ggplot_req_for <- function(pkg) {
  seen <- character()
  queue <- c(pkg, all_deps(pkg))
  queue <- unique(queue[queue %in% rownames(ap)])
  
  reqs <- list()  # keep all hits with source package
  
  for (p in queue) {
    fields <- dep_strings(p)
    deps_df <- parse_deps(fields)
    if (!nrow(deps_df)) next
    hit <- subset(deps_df, tolower(pkg) != "ggplot2" & tolower(pkg) == tolower(pkg)) # placeholder no-op
    # filter for ggplot2 rows
    gg <- subset(deps_df, tolower(pkg) == "ggplot2")
    if (nrow(gg)) {
      reqs[[p]] <- gg$constraint
    }
  }
  
  if (!length(reqs)) return(list(max_needed = NA_character_, sources = character()))
  # normalize constraints like ">= 3.5.0"
  nums <- sapply(reqs, function(v) {
    v <- v[nzchar(v)]
    if (!length(v)) return(NA_character_)
    m <- regmatches(v, regexpr(">=\\s*([0-9.]+)", v))
    if (!length(m)) return(NA_character_)
    sub(">=\\s*", "", m)
  })
  
  nums <- unlist(nums, use.names = TRUE)
  nums <- nums[!is.na(nums)]
  
  if (!length(nums)) return(list(max_needed = NA_character_, sources = names(reqs)))
  
  # strongest requirement = max version among constraints
  ord <- order(package_version(nums), decreasing = TRUE)
  list(max_needed = nums[ord][1], sources = names(nums)[ord])
}

cat("\n=== Checking suspects against ggplot2", target, "===\n")
culprits <- list()

for (s in suspects) {
  if (!s %in% rownames(ap)) {
    cat(sprintf("… %-10s: not found on CRAN index\n", s))
    next
  }
  info <- ggplot_req_for(s)
  req  <- info$max_needed
  if (is.na(req)) {
    cat(sprintf("✓ %-10s: no explicit ggplot2 constraint found (direct or deps)\n", s))
  } else {
    ok <- utils::compareVersion(req, target) <= 0
    mark <- if (ok) "OK" else "FORCES UPGRADE"
    cat(sprintf("%s %-10s: requires ggplot2 >= %s\n", if (ok) "✓" else "✗", s, req))
    if (!ok) {
      culprits[[s]] <- list(required = req, sources = info$sources)
      # print first couple of sources
      srcs <- paste(head(info$sources, 3), collapse = ", ")
      cat("    via:", srcs, if (length(info$sources) > 3) ", …", "\n")
    }
  }
}

if (length(culprits)) {
  cat("\n=== Culprit summary (need pins/downgrade) ===\n")
  for (nm in names(culprits)) {
    cat(sprintf("- %s → ggplot2 >= %s\n", nm, culprits[[nm]]$required))
  }
} else {
  cat("\nNo suspect explicitly requires ggplot2 > ", target, ".\n",
      "Then the upgrade is likely coming from a NON-suspect dependency (e.g., tidyverse branch) or a stale lockfile.\n", sep = "")
}
