# Topographic structuring of arthropod communities revealed by leaf-derived environmental DNA on Oʻahu, Hawaiʻi

Analysis code and processed ASV data for the accompanying study: arthropod communities along five Oʻahu ridges, with comparisons among ʻōhiʻa (*Metrosideros polymorpha*), koa (*Acacia koa*) and strawberry guava (*Psidium cattleyanum*) at ʻAiea and Mānana. Native/introduced classifications describe biogeographic origin in the Hawaiian archipelago.

**Raw FASTA sequence data:** download `raw_files.zip` from [Dryad (10.5061/dryad.2z34tmq09)](https://doi.org/10.5061/dryad.2z34tmq09).

## Analysis files

| File in `scripts/` | Analysis and manuscript outputs |
| --- | --- |
| [00_source.R](scripts/00_source.R) | Shared inputs, control adjustment, sample subsets and plotting helpers. |
| [01_map.qmd](scripts/01_map.qmd) | Sampling map (Figure 1). |
| [02_initial_analyses.qmd](scripts/02_initial_analyses.qmd) | NMDS, PERMANOVA, dispersion, turnover/nestedness and distance–decay (Figures S3–S4). |
| [03_overlap_invasion.qmd](scripts/03_overlap_invasion.qmd) | ASV overlap and introduced-ASV proportions (Figures 4–5; S2, S8–S9). |
| [04_species_richness.qmd](scripts/04_species_richness.qmd) | Richness, elevation and rainfall (Figures S6–S7). |
| [05_titan_breakpoints.qmd](scripts/05_titan_breakpoints.qmd) | TITAN2 and presence/absence comparison (Figure 6; S10). |
| [06_molecular_framing.qmd](scripts/06_molecular_framing.qmd) | Reference-based and NIClassify status summaries (Figure 2; Table 1; S1). |
| [07_status_communities.qmd](scripts/07_status_communities.qmd) | Native/introduced community structure and within-ridge distance–decay (Figure 3; S5). |
| [server/ohia_gradient.sh](scripts/server/ohia_gradient.sh) | ANML sequence-processing command record. |
| [server/trimflow.py](scripts/server/trimflow.py) | Primer/barcode trimming with the bundled Bash engine and configuration. |
| [server/blast2taxonomy_v1.4.2.py](scripts/server/blast2taxonomy_v1.4.2.py) | Taxonomy annotation by [Yannis Schöneberg](https://github.com/yschoeneberg/blast2taxonomy) (full documentation); [CC0 license](scripts/server/blast2taxonomy_LICENSE). |

## Inputs

- [ohia_gradient_2026-02-09.xlsx](data/ohia_gradient_2026-02-09.xlsx): ASV read tables, sample metadata and taxonomy; analyses read `zotu_table_merged`, `metadata` and `taxonomy`.
- [NIClassify predictions](NIClassify/06_NIClassify_predict.csv): saved classifier results and reference-status information.
- [tax_invasive.xlsx](data/tax_invasive.xlsx): curated reference assignments (`taxonomy` and `taxonomy_bold`).
- [Elevation raster](data/oahu_elevation_clipped.tif) and [tree illustrations](data/images/): figure backgrounds and icons.
- [ANML primer table](scripts/server/primer_database.csv): the paper's primer sequences and expected insert length.

## Run the R analyses

Use **R 4.5.1**, **Quarto** and Chrome/Chromium for table and Sankey exports. Restore the recorded packages, including **ggplot2 3.5.2**, from the repository root:

```sh
Rscript -e 'renv::restore()'
make figures
```

Run individual notebooks with `quarto render scripts/03_overlap_invasion.qmd`.

Main figures and Table 1 go to `figures/main/`; supplementary figures go to `figures/supp/`; temporary Sankey files go to `output/`.

## Server sequence processing

Create the trimming environment from the repository root:

```sh
mamba env create -f scripts/server/environment.yaml
mamba activate trimming_env
```

Launch [trimflow.py](scripts/server/trimflow.py) from the directory containing merged `.fasta` reads. Select **ANML**, **A,B,C**, **UNIX** and the output directory. It uses the bundled ANML primer table by default. The [command record](scripts/server/ohia_gradient.sh) lists the other processing tools, original reads and reference database required for the complete upstream workflow.
