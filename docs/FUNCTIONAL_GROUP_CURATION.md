# Functional-group curation

## Purpose and scope

The 13 functional groups were predefined to test hypothesized direct sunlight responses and substrate-mediated legacies relevant to litter decomposition. They are a targeted subset of microbial functional potential, not a universal metabolic ontology.

## Treatment-blind candidate review

The candidate universe was the union of all 7,931 KOs detected at nonzero abundance in either final experiment. Candidate discovery used four complementary searches of KO names and structured KEGG pathway, module, BRITE, and entry information. Searches were intentionally sensitive and identified records for review; keyword or pathway membership never assigned a KO automatically.

Candidate assessment considered the current documented molecular function, substrate specificity, direct versus auxiliary involvement, and the prespecified group boundary. KOs without direct molecular-function evidence were excluded even when broad pathway co-membership produced a search hit. The full universe table preserves the screen hits, evidence snapshot, disposition, and rationale for both inclusions and exclusions.

Candidate evidence and draft rationales were organized computationally. Final inclusion, exclusion, and group-assignment decisions were reviewed by the authors. Treatment identities and statistical results were not used to decide membership.

## Assignment rule

Each of the 332 retained KOs has one primary functional-group assignment. This prevents double counting in group abundance sums but simplifies multifunctional biology. The complete mapping and evidence are distributed as `data/functional_metagenomics/frozen_332_KO_functional_group_mapping.csv` and Supplementary Table S2.

## Sensitivity and secondary analysis

Primary group-level results are accompanied by log-transformation and leave-one-sample-out sensitivity analyses. An experiment-specific 75%-prevalence KO-level CLR screen is retained as exploratory evidence; it does not redefine the prespecified groups and does not supersede the group-level endpoint.
