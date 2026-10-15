.PHONY: figures

figures:
	quarto render scripts/01_map.qmd
	quarto render scripts/02_initial_analyses.qmd
	quarto render scripts/03_overlap_invasion.qmd
	quarto render scripts/04_species_richness.qmd
	quarto render scripts/05_titan_breakpoints.qmd
	quarto render scripts/06_molecular_framing.qmd
	quarto render scripts/07_status_communities.qmd
