# Software versions

Release test platform:

- Linux (x86_64)
- R 4.1.2
- Python 3.13.2

Analysis package versions are pinned in `R-packages.tsv` and `requirements.txt`. The upstream, nonredistributed bioinformatics workflow used BBDuk/BBMap, Bowtie2, SAMtools, HUMAnN 3.9, and DIAMOND; exact upstream command parameters belong in the manuscript Methods. Database provenance is recorded separately in `docs/DATABASE_VERSIONS.md`.

The package was run successfully both from its project location and from an isolated copied directory. Each R workflow writes a session-information file into `results/reproduced/`.
