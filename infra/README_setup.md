Use renv::restore() first

If ggplot breaks (4.x), run one of:

source("scripts/00_renv_ggplot35x.R") (preferred)

or the “force” script (as fallback)

This removes tribal knowledge and prevents future-you from re-debugging ggplot at 2 a.m.
