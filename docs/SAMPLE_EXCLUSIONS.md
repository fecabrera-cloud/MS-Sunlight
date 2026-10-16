# Sample inclusion and exclusions

The authoritative machine-readable record is `data/sample_inclusion_ledger.csv`.

Key principles:

- exclusions are assay specific;
- the confirmed field extraction failure is removed consistently from shotgun, 16S, and ITS representations before downstream analysis;
- marker libraries below 15,000 target-marker reads after taxonomy filtering are excluded before prevalence filtering and rarefaction;
- differences in final sample counts among assays are retained rather than forcing artificial intersection across unrelated endpoints;
- missing measurements remain missing and are not imputed.

Final primary sample sizes reproduced by the workflow are:

| Analysis | Field | Microcosm |
|---|---:|---:|
| 16S amplicons | 19 | 36 |
| ITS amplicons | 18 | 40 |
| Shotgun functional metagenomics | 19 | 41 |
| Functional measurements | — | 46 for dry-mass loss, biomass, and beta-glucosidase; 43 for biomass/respiration ratio |
