# Figure and table map

| Item | Final file under figures/final | Reproduced output under results/reproduced | Generating script under scripts |
|---|---|---|---|
| Figure 1 | `Figure1.pptx` | `Author artwork` | `—` |
| Figure 2 | `Figure2.pdf` | `amplicons/field_amplicon_NMDS.pdf` | `run_amplicon_audit.R` |
| Figure 3 | `Figure3.pdf` | `functional_metagenomics/field/pcoa_publication.pdf` | `run_functional_metagenomics_audit.R` |
| Figure 4 | `Figure4.pdf` | `amplicons/Figure4.pdf` | `run_amplicon_audit.R` |
| Figure 5 | `Figure5.pdf` | `functional_metagenomics/microcosm/pcoa_publication.pdf` | `run_functional_metagenomics_audit.R` |
| Figure 6 | `Figure6.pdf` | `functional_measurements/Figure6.pdf` | `run_functional_measurements_audit.R` |
| Figure S1 | `FigS1.pdf` | `microclimate/temperature_day_night_boxplots.pdf` | `run_microclimate_analysis.py` |
| Figure S2 | `FigS2.pdf` | `microclimate/temperature_day_night_timeseries.pdf` | `run_microclimate_analysis.py` |
| Figure S3 | `FigS3.pdf` | `microclimate/humidity_day_night_boxplots.pdf` | `run_microclimate_analysis.py` |
| Figure S4 | `FigS4.pdf` | `microclimate/humidity_day_night_timeseries.pdf` | `run_microclimate_analysis.py` |
| Figure S5 | `FigS5.pdf` | `amplicons/field_top12_genera.pdf` | `run_amplicon_audit.R` |
| Figure S6 | `FigS6.pdf` | `amplicons/field_all_genus_volcano.pdf` | `run_amplicon_deseq.py` |
| Figure S7 | `FigS7.pdf` | `functional_metagenomics/field/supplementary_functional_groups.pdf` | `run_functional_metagenomics_audit.R` |
| Figure S8 | `FigS8.pdf` | `amplicons/microcosm_top12_genera.pdf` | `run_amplicon_audit.R` |
| Figure S9 | `FigS9.pdf` | `amplicon_comparison/field_microcosm_genus_comparison.pdf` | `plot_genus_comparison.R` |
| Figure S10 | `FigS10.pdf` | `amplicons/microcosm_all_genus_volcano.pdf` | `run_amplicon_deseq.py` |
| Figure S11 | `FigS11.pdf` | `functional_metagenomics/microcosm/supplementary_functional_groups.pdf` | `run_functional_metagenomics_audit.R` |

S8 retains one bar per sample, ordered by substrate and inoculum. S9 is descriptive. S10 retains the existing genus-level analysis. S11 contains only the functional-group overview.

Tables S1–S2 are produced by `build_supplementary_tables.py`. The figure bundle is produced by `build_figure_bundle.py` from the final PDFs. The final PDF/PPTX assets are protected from overwrite during reproduction.
