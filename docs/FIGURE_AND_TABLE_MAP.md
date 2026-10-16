# Figure and table map

| Manuscript item | Final asset | Reproduced output | Generating script |
|---|---|---|---|
| Figure 1 | `figures/final/Figure1.pptx` | Author artwork | — |
| Figure 2 | `figures/final/Figure2.pdf` | `results/reproduced/amplicons/field_amplicon_NMDS.*` | `scripts/run_amplicon_audit.R` |
| Figure 3 | `figures/final/Figure3.pdf` | `results/reproduced/functional_metagenomics/field/pcoa_publication.*` | `scripts/run_functional_metagenomics_audit.R` |
| Figure 4 | `figures/final/Figure4.pdf` | `results/reproduced/amplicons/microcosm_amplicon_NMDS.*` | `scripts/run_amplicon_audit.R` |
| Figure 5 | `figures/final/Figure5.pdf` | `results/reproduced/functional_metagenomics/microcosm/pcoa_publication.*` | `scripts/run_functional_metagenomics_audit.R` |
| Figure 6 | `figures/final/Figure6.pdf` | `results/reproduced/functional_measurements/primary_functional_responses.*` | `scripts/run_functional_measurements_audit.R` |
| Figure S1 | `figures/final/FigS1.pdf` | `results/reproduced/amplicons/field_top12_genera.*` | `scripts/run_amplicon_audit.R` |
| Figure S2 | `figures/final/FigS2.pdf` | `results/reproduced/amplicons/field_all_genus_volcano.*` | `scripts/run_amplicon_deseq.py` |
| Figure S3 | `figures/final/FigS3.pdf` | `results/reproduced/functional_metagenomics/field/supplementary_functional_groups.*` | `scripts/run_functional_metagenomics_audit.R` |
| Figure S4 | `figures/final/FigS4.pdf` | `results/reproduced/amplicons/microcosm_top12_genera.*` | `scripts/run_amplicon_audit.R` |
| Figure S5 | `figures/final/FigS5.pdf` | `results/reproduced/amplicons/microcosm_all_genus_volcano.*` | `scripts/run_amplicon_deseq.py` |
| Figure S6 | `figures/final/FigS6.pdf` | `results/reproduced/functional_metagenomics/microcosm/supplementary_functional_groups.*` | `scripts/run_functional_metagenomics_audit.R` |
| Figure S7 | `figures/final/FigS7.pdf` | `results/reproduced/microclimate/temperature_day_night_boxplots.*` | `scripts/run_microclimate_analysis.py` |
| Figure S8 | `figures/final/FigS8.pdf` | `results/reproduced/microclimate/temperature_day_night_timeseries.*` | `scripts/run_microclimate_analysis.py` |
| Figure S9 | `figures/final/FigS9.pdf` | `results/reproduced/microclimate/humidity_day_night_boxplots.*` | `scripts/run_microclimate_analysis.py` |
| Figure S10 | `figures/final/FigS10.pdf` | `results/reproduced/microclimate/humidity_day_night_timeseries.*` | `scripts/run_microclimate_analysis.py` |
| Table 1 | Manuscript | `results/reproduced/chemistry/Table_1_descriptive_no_inference.csv` | `scripts/run_chemistry_analysis.py` |
| Table S1 | `tables/Table_S1_functional_group_definitions.csv` | Same | `scripts/build_supplementary_tables.py` |
| Table S2 | CSV, XLSX, and DOCX under `tables/` | Final KO mapping | `scripts/build_supplementary_tables.py` |
| Table S3 | CSV, XLSX, and DOCX under `tables/` | Functional sensitivity output | `scripts/build_supplementary_tables.py` |

Running `scripts/run_all.sh` creates the `results/` and `validation/` directories locally. Both directories are ignored by Git and are not part of the public repository.
