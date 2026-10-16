# Analysis boundary

## Included

This repository begins with processed feature tables:

- taxonomically annotated 16S rRNA and ITS ASV count tables;
- unstratified KO relative-abundance matrices produced by HUMAnN regrouping;
- the frozen 7,931-KO review universe and final 332-KO/13-group mapping;
- cleaned microcosm functional measurements;
- replicate-level substrate-chemistry measurements;
- anomaly-filtered field microclimate observations.

All manuscript-level statistics and base visualizations from these inputs are reproduced.

## Not included

- raw FASTQ data;
- adapter/quality trimming and host-removal workflows;
- QIIME 2 denoising/taxonomic-classification artifacts;
- per-sample HUMAnN execution and UniRef90 database files;
- Slurm/HPC orchestration;
- third-party KEGG or UniProt database distributions.

Those upstream steps are described in the manuscript. Raw-data archive identifiers belong in `RAW_DATA_ACCESSIONS.md`.

## Statistical units and constraints

- Field community and functional analyses use field replicate/block to restrict permutations or as a fixed blocking term, as appropriate.
- Microcosm community and multivariate functional analyses restrict permutations within inoculum-by-replicate blocks.
- Microcosm univariate functional-group and response models use inoculum as a random intercept and account for replicate/block as specified by each script.
- Chemistry values are descriptive because repeat readings are analytical, not independent treatment replicates.
- The biomass-to-respiration ratio is reported descriptively as calculated and is not treated as a direct physiological measurement of microbial carbon-use efficiency.

## Final figures versus reproduced plots

`results/reproduced/` contains plots generated directly by the analysis code. `figures/final/` contains manuscript-ready author assets. This separation is deliberate: minor Illustrator or layout edits are not falsely presented as programmatically reproduced, and rerunning the statistics cannot overwrite the publication figures.
