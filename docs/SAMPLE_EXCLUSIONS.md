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
| Functional measurements | — | 43 for dry-mass loss, biomass, and beta-glucosidase; 40 for biomass/respiration ratio |

Functional-measurement analyses are restricted to the five experimental blocks. Three additional full-sunlight bottles (IDs 46–48, one per inoculum origin, coded block 6) were incubated alongside the experiment but were not processed for molecular analyses. They remain in the cleaned measurement dataset with `in_primary_analysis = FALSE` and are excluded from the fitted models, so that functional measurements and molecular analyses share the same five-block design. The retained inoculated records comprise 15 full-sunlight, 15 reduced-UV, and 13 reduced-UV-and-blue-green substrates. Bottle IDs 29 and 35 have no record in this dataset, so the analysed sample size is 43 rather than the 45 bottles initially established. The seven non-inoculated controls remain available for descriptive purposes and are excluded from the fitted models. Cumulative respiration is unavailable for retained bottle IDs 9, 44, and 45; the biomass-to-respiration ratio therefore has N = 40. No values were imputed.
