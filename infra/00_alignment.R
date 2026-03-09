# Load necessary libraries
libs <- c("tidyverse","openxlsx","janitor", "vegan", "ggtext", "RColorBrewer",
          "glue", "here", "ggimage", "scales", "viridis","Biostrings", "ape", "phangorn", "parallel",
          "ggtree", "ape", "msa", "seqinr", "seqinr", "DECIPHER", "phangorn")
# Install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}

# Load libraries
invisible(lapply(libs, library, character.only = T))

rm(installed_libs)


tax_unfiltered <- openxlsx::read.xlsx(here::here("data", "ANML.xlsx"), sheet = "tax_unfiltered") |> 
  janitor::clean_names() |>
  filter(phylum == "Arthropoda") |>
  filter(sbjct_len > 170)


tax_filtered <- tax_unfiltered |> 
  group_by(query, kingdom, phylum, order, family, genus, species) |>
  count() |> 
  ungroup() |>
  group_by(query) |>
  slice_max(n = 1, order_by = n, with_ties = FALSE) |> 
  ungroup() |>
  mutate(otu_num = as.numeric(str_replace(query, "Zotu", "")),.before =2) |>
  arrange(otu_num) |> 
  as.data.frame()


sequences_df_filtered <- sequences_df |>
  inner_join(tax_filtered, by = c("zotu_id" = "query"))


#create alignment with MSA based on clustalW
alignment_zotus <- msa::msa(sequences_df_filtered $sequence, type = "dna", order = "input", method = "ClustalW")

# Convert alignment to a DNAbin object
alignment_dnabin <- as.DNAbin(alignment_zotus)

# Convert alignment to phyDat for filtering
alignment_phy <- as.phyDat(alignment_dnabin)

# Remove sequences with too many gaps
filtered_alignment <- alignment_phy[sapply(alignment_phy, function(x) sum(x != "-") > 10)] 

# Convert back to DNAbin
filtered_dnabin <- as.DNAbin(alignment_zotus)

# Save alignment to a file (e.g., PHYLIP format)
write.dna(filtered_dnabin, file = here::here("output/aligned_zotus.phy"), format = "sequential")
write.dna(alignment_dnabin, file = here::here("output/aligned_zotus.fas"), , format = "fasta")

# Compute distance matrix
dist_matrix <- dist.dna(filtered_dnabin, model = "TN93")  # Tamura-Nei model

# Create new tip labels by combining zOTU IDs with species names
custom_labels <- paste0(seqs$zotu_id, "_", seqs$species)

# Convert DNAbin matrix to list format (preserves labels correctly)
filtered_dnabin <- as.list(filtered_dnabin)

# Assign new names to sequences
names(filtered_dnabin) <- custom_labels


# Build NJ tree
nj_tree <- njs(dist_matrix)

# Update tip labels in the tree
nj_tree$tip.label <- custom_labels

# Plot NJ tree
plot(nj_tree, main = "Neighbor-Joining Tree of zOTUs")


# Build ML tree
tree_start <- nj_tree  # If this still fails, skip NJ
fit <- pml(tree_start, data = as.phyDat(filtered_dnabin))
fit_opt <- optim.pml(fit, model = "GTR", optGamma = TRUE, optInv = TRUE)
plot(fit_opt$tree, main = "ML Tree")

# Set number of bootstrap replicates
bs_replicates <- 100  # Increase to 1000 for publication-quality results

# Run bootstrap analysis
bs_tree <- bootstrap.pml(fit_opt, bs = bs_replicates, optNni = TRUE)

# Save bootstrap values
tree_bootstrapped <- fit_opt$tree
tree_bootstrapped$node.label <- bs_tree

# Save tree to file (optional)
write.tree(tree_bootstrapped, file = here("output/ML_tree_bootstrapped.nwk"))

# Plot ML tree with bootstrap values and order nodes
ggtree(tree_bootstrapped, layout = "rectangular") + 
  geom_tiplab(size = 2, align = TRUE)   # Label taxa

# Ensure we have the correct zOTU-to-Order mapping
zotu_order_map <- setNames(seqs$order, seqs$zotu_id)

# Convert tree to data frame
bs_tree_data <- fortify(tree_bootstrapped)

# Convert tree to data frame
bs_tree_data <- fortify(tree_bootstrapped) |> 
  mutate(zotu_id_extracted = ifelse(isTip, str_extract(label, "Zotu[0-9]+"), NA))  # Extract zOTU IDs only for tips

bs_tree_data_combined <- bs_tree_data %>%
  mutate(order = zotu_order_map[zotu_id_extracted]) 



ggtree(bs_tree_data_combined, layout = "rectangular", aes(color = order)) +
  geom_tiplab(aes(label = label), size = 2, align = TRUE) +  # Keep zOTU_Species labels
  #scale_color_manual(values = order_colors, na.value = "black") +  # Custom colors for orders
  theme_minimal() +
  ggtitle("ML Tree with Branch Colors by Order")




