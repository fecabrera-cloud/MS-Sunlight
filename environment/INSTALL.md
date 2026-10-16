# Environment installation

The release was validated with R 4.1.2 and Python 3.13.2 on Linux. Exact package versions are recorded in `R-packages.tsv` and `requirements.txt`.

Python packages can be installed in a clean virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/requirements.txt
```

R packages are available from CRAN except `phyloseq`, which is distributed through Bioconductor. To install current compatible versions:

```r
install.packages(c(
  "dplyr", "emmeans", "ggplot2", "ggrepel", "multcompView", "nlme",
  "patchwork", "permute", "scales", "stringr", "tibble", "tidyr", "vegan"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
```

Because archived binary availability depends on operating system and R version, exact reconstruction is best performed with a container or package snapshot built from the recorded versions. Each R workflow also writes `R_session_info.txt` beside its results.
