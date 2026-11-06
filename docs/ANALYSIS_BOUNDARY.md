# Analysis boundary

The repository begins with processed ASV counts and taxonomy, unstratified HUMAnN KO profiles, the frozen KO mapping, cleaned functional measurements, substrate chemistry, and field microclimate observations. Raw FASTQ processing, denoising, host removal, HUMAnN execution, HPC orchestration, and third-party databases are outside its scope.

See [ANALYSIS_METHODS.md](ANALYSIS_METHODS.md) for the model specifications, experimental units, and interpretation limits. Chemistry is summarized descriptively because analytical measurements do not independently replicate field treatment pools.

Analysis scripts write regenerated data, models, diagnostics, and plots under `results/reproduced/`. Final PDF/PPTX assets under `figures/final/` are protected from changes by the end-to-end runner. The supplementary Word figure bundle is derived from those final PDFs and may be rebuilt. The manuscript document is maintained separately.
