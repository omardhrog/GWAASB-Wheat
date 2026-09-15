# GWAASB-Wheat

A reproducible **grain-yield(`GY`)** workflow for integrating genomic prediction and stability analysis in multi-environment wheat trials.

## Scope

The pipeline implements:

- environment-specific BLUE estimation;
- genomic relationship matrix construction;
- GBLUP and GBLUP-G×E models;
- leakage-controlled CV1 and CV2;
- within-environment prediction accuracy;
- genomic G×E-effect matrix extraction;
- GWAASB calculation by singular-value decomposition;
- comparison with phenotypic WAASB, GE row standard deviation, and Wricke's ecovalence;
- joint assessment of predicted genomic performance and stability.

## Privacy boundary

This repository is intentionally restricted to **grain yield (`GY`)**.

## Required input files

Place these files in the R working directory (normally the repository root), or edit `POP_FILES`:

| File | Required structure |
| --- | --- |
| `phenotype_data_DHLS-1.csv` | `ENV, GEN, REP, GY` |
| `phenotype_data_DHLS-2.csv` | `ENV, GEN, REP, GY` |
| `SNP_clean_fixed_DHLS-1.csv` | `GEN` followed by SNP-marker columns |
| `SNP_clean_fixed_DHLS-2.csv` | `GEN` followed by SNP-marker columns |

SNP genotypes may be coded as `0/1/2` or `-1/0/1`. Missing marker values are mean-imputed when the genomic relationship matrix is constructed.

## Run

```r
source("R/GWAASB_Wheat_GY_pipeline.R")
```

The default analysis uses two populations, excludes parents `G101` and `G102`, retains 100 DH lines per population, and performs 100 repetitions of five-fold CV1 and CV2. These settings can be changed in the **SETTINGS** section.

## Outputs

All generated files are written to `results/`, including marker-QC summaries, environment BLUEs, genomic predictions, GWAASB values, stability comparisons, CV summaries, environment-specific metrics, validation predictions, and session information.

## Reproducibility note

Cross-validation partitions are created before phenotype adjustment. Validation observations are adjusted using training-derived replicate effects, and the same partitions are used for GBLUP and GBLUP-G×E within each repetition.

## Status

Research code associated with the GWAASB wheat analysis. The repository is currently private and contains no study data.
