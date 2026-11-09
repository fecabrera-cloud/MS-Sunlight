# Data dictionary

## Amplicons

### `data/amplicons/16S_feature_table.csv` and `ITS_feature_table.tsv`

Rows are ASVs and columns are paper-scope samples. Values are nonnegative read counts after upstream denoising and taxonomic processing. The first column is the ASV identifier.

### Sample metadata

- `Experiment` / `experiment`: field inoculum endpoint (`Inoculos`) or reciprocal microcosm (`RecyprocalT`).
- `TreatmentA` / `treatA`: field treatment or microcosm substrate photo-history.
- `TreatmentB` / `treatB`: field litter identity or microcosm inoculum origin.
- `Replicate` / `replicate`: experimental replicate/block.
- `sample_id` / `sample-id`: marker-specific library identifier.

### Taxonomy

Taxonomy strings contain semicolon-separated ranks. The scripts remove non-target and organellar ASVs, apply a 5% within-analysis prevalence filter, and rarefy only for Bray-Curtis ordination/PERMANOVA. Differential abundance uses unrarefied genus counts.

## Functional metagenomics

### KO matrices

`field_KO_relative_abundance.tsv` and `microcosm_KO_relative_abundance.tsv` contain samples in rows and unstratified KOs in columns. Values are total-sum-scaled relative abundances derived from HUMAnN RPK profiles. The denominator is the complete retained KO profile, not only the targeted subset.

### `frozen_332_KO_functional_group_mapping.csv`

One row per retained KO. Key fields include:

- `ko_id`, `kegg_name`: KO identifier and functional name;
- `functional_group`: sole primary analysis group;
- `hypothesis_role`: direct sunlight response or substrate-mediated role;
- `assignment_confidence`, `curation_reason`, `evidence_basis`: assignment evidence and rationale;
- pathway/module/BRITE evidence and KEGG entry URL;
- experiment-specific detection, prevalence, and mean relative abundance;
- `curation_note`: additional context for boundary-sensitive assignments.

### `functional_group_definitions.csv`

One row per targeted group with scope, inclusion/exclusion boundaries, hypothesis role, and mapped/detected KO counts.

### `full_7931_KO_universe_review.csv`

Treatment-blind audit of every KO detected at nonzero abundance in either final experiment. It records four candidate-discovery screens, KEGG evidence snapshots, final inclusion/exclusion decisions, and rationales. An exclusion denotes lack of direct molecular-function evidence for the predefined groups, not lack of biological importance.

## Functional measurements

`cleaned_functional_measurements.csv` contains inoculated microcosms and noninoculated controls.

The dataset contains 43 inoculated bottles from experimental blocks 1–5 and seven non-inoculated controls, coded as block 7. Records labelled block 6 were excluded; see `docs/SAMPLE_EXCLUSIONS.md`. Three retained inoculated bottles lack cumulative respiration, leaving 40 biomass-to-respiration ratios.

Primary variables are:

- `decomposition_pct`: litter dry-mass loss during incubation; microbial biomass was not subtracted;
- `microbial_biomass_c_ug_g`: biomass C inferred from substrate-induced respiration;
- `beta_glucosidase`: potential beta-glucosidase activity;
- `respiration_auc`: integrated respiration response;
- `cue_proxy`: microbial biomass C divided by integrated respiration; three values are missing because cumulative respiration was unavailable;
- `inoculum`, `block`, `substrate_code`: model design variables.

## Chemistry

`combined_replicate_level_source_data.csv` is a long-format integration of chemical, fiber, and muffle-furnace measurements. `assay_unit_id` identifies the analytical unit within each source assay. Chemistry summaries report assay-specific sample sizes, means, and standard deviations by substrate treatment.

## Microclimate

`field_only_anomaly_filtered.csv` contains logger timestamp, temperature, humidity, treatment, logger replicate/serial, and day/night period after declared anomaly filtering. Removed values and reasons are preserved separately.

## Sample ledger

`data/sample_inclusion_ledger.csv` consolidates assay-specific inclusion/exclusion information. Sample counts legitimately differ across markers and measurements because technical failures and sequencing thresholds were assay specific.
