# Topographic structuring of arthropod communities revealed by leaf-derived environmental DNA on Oʻahu, Hawaiʻi

Analysis code and processed inputs for the accompanying study of arthropod communities detected from foliage eDNA. The study compares communities along five Oʻahu ridges and among ʻōhiʻa (*Metrosideros polymorpha*), koa (*Acacia koa*) and strawberry guava (*Psidium cattleyanum*) at ʻAiea and Mānana.

The R workflow starts from processed ASV tables and saved NIClassify predictions. Native and introduced refer to biogeographic origin in the Hawaiian archipelago; an introduced classification does not establish invasiveness. Sequence detections are not direct measures of organism abundance.

## Analysis files

| File in `scripts/` | Analysis and manuscript outputs |
| --- | --- |
| [00_source.R](scripts/00_source.R) | Shared data import, control adjustment, sample subsets and plotting helpers. |
| [01_map.qmd](scripts/01_map.qmd) | Sampling map (Figure 1). |
| [02_initial_analyses.qmd](scripts/02_initial_analyses.qmd) | NMDS, PERMANOVA, dispersion, turnover/nestedness and spatial distance–decay (Figures S3–S4). |
| [03_overlap_invasion.qmd](scripts/03_overlap_invasion.qmd) | ASV overlap, status comparisons and introduced-ASV proportions (Figures 4–5; S2, S8–S9). |
| [04_species_richness.qmd](scripts/04_species_richness.qmd) | Richness, elevation and rainfall relationships (Figures S6–S7). |
| [05_titan_breakpoints.qmd](scripts/05_titan_breakpoints.qmd) | TITAN2 elevation analysis and presence/absence sensitivity comparison (Figure 6; S10). |
| [06_molecular_framing.qmd](scripts/06_molecular_framing.qmd) | Reference-based and NIClassify status summaries (Figure 2; Table 1; S1). |
| [07_status_communities.qmd](scripts/07_status_communities.qmd) | Native/introduced community structure and within-ridge distance–decay (Figure 3; S5). |
| [ohia_gradient.sh](scripts/ohia_gradient.sh) | Commented record of the upstream ANML sequence-processing commands. |

## Inputs

- `data/ohia_gradient_2026-02-09.xlsx`: sample metadata, ASV read tables and taxonomy; shared analyses read the `metadata`, `zotu_table_merged` and `taxonomy` sheets.
- `NIClassify/06_NIClassify_predict.csv`: saved classifier predictions and reference-status information.
- `data/tax_invasive.xlsx`: curated reference-based status assignments (`taxonomy` and `taxonomy_bold` sheets).
- `data/oahu_elevation_clipped.tif` and `data/images/`: elevation background and tree illustrations used in the figures.

## Run the R analyses

Use **R 4.5.1**, **Quarto** and the package versions recorded in `renv.lock`, including **ggplot2 3.5.2**. Chrome or Chromium is needed for table and Sankey image exports. Run from the repository root:

```sh
Rscript -e 'renv::restore()'
make figures
```

To run one analysis:

```sh
quarto render scripts/03_overlap_invasion.qmd
```

Figures are written to `figures/main/` and `figures/supp/`; derived tables and temporary graphics go to `output/`. Filenames include the execution date. Permutation, bootstrap and ordination results can vary where the original code does not fix the random seed.

The upstream shell record requires the original reads, external tools, database and interactive `trimflow.py` implementation. It documents the recovered ANML commands and is separate from the processed-data R workflow; see its comments for requirements.
