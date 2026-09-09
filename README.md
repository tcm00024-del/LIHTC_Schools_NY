# Low Income Housing and School Performance: Evidence from New York State Report Cards

Replication data and code for the paper examining the effects of Low-Income Housing Tax Credit (LIHTC) funded developments on nearby New York State public schools.

**Author:** Thaddeus Meadows, Department of Economics, West Virginia University
**Contact:** tcm00024@mix.wvu.edu

## Overview

This repository contains the analysis-ready data, data-cleaning pipeline, and estimation code supporting the paper. The paper investigates whether the opening of LIHTC-funded housing developments affects standardized test scores, student body composition, and attendance/suspension outcomes at nearby New York State public schools (excluding New York City), using Callaway-Sant'Anna and Sun-Abraham staggered-treatment estimators.

## Repository structure

```
data/
  processed/    Analysis-ready panel data used directly by the R scripts in code/estimation
  README.md     Data dictionary and description of raw source data (not redistributed here)
code/
  cleaning/     Python/Google Colab notebooks that build the processed panel from raw HUD and NYSED source files
  estimation/   R scripts implementing TWFE, Callaway-Sant'Anna, and Sun-Abraham event study estimators
    enrollment/ R scripts for the enrollment/demographic-composition outcomes (headcount and
                demographic panels); shares helper functions via enrollment_estimation_helpers.R
    scores/     R scripts for the test-score outcomes, by subgroup (Full and No-NYC samples,
                each split into Full/ and No_NYC/ subfolders under data/processed/Subgroups/);
                shares helper functions via cs_estimation_helpers.R
```

## Data sources

- **HUD LIHTC administrative data**: U.S. Department of Housing and Urban Development, Office of Policy Development and Research. Publicly available at https://lihtc.huduser.gov/
- **NYSED school-level data**: New York State Education Department school report card data (test scores, enrollment, attendance, suspensions). Publicly available at https://data.nysed.gov/
- **School district shapefiles**: publicly available ArcGIS shapefiles for New York public school districts.

Raw source files are not redistributed in this repository; see `data/README.md` for details on obtaining them and for a full data dictionary of the processed panel.

## Requirements

- Python 3.x (Google Colab environment) for the data-cleaning pipeline
- R with the following packages for estimation: `did` (Callaway-Sant'Anna), `fixest` (Sun-Abraham via `sunab()` and TWFE via `feols()`), `data.table`, `ggplot2` (event-study plots), `here` (portable, project-root-relative file paths — see note below)

All estimation scripts locate the repository root automatically via the `here` package rather than a hardcoded path. This works out of the box after `git clone`; if you instead download the repository as a ZIP (no `.git` folder), add an empty file named `.here` at the repository root first so `here()` can still find it.

## Reproducing the analysis

1. Run the notebooks in `code/cleaning/` in numbered order to reconstruct the processed panel from raw source files (see `data/README.md` for raw data access).
2. Run the scripts in `code/estimation/` against the processed panel in `data/processed/` to reproduce the paper's tables and figures.

## Citation

If you use this data or code, please cite the paper and/or this repository. See `CITATION.cff` for citation metadata, or use the "Cite this repository" link on GitHub.

## License

Code in this repository is licensed under the MIT License (see `LICENSE`). Data files in `data/processed/` are licensed under CC-BY 4.0.
