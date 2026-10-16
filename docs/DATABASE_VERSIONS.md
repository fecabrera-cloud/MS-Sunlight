# Database and annotation provenance

## Confirmed

- Functional profiling software: HUMAnN 3.9.
- Translated search: DIAMOND against UniRef90.
- UniRef90 archive recorded in project provenance: `uniref90_annotated_v201901b_full.tar.gz`.
- Regrouping command family: `humann_regroup_table --groups uniref90_ko`.
- Utility-mapping filename: `map_ko_uniref90.txt.gz`.
- Downstream analysis uses only unstratified KO rows; `UNMAPPED` and `UNGROUPED` are not biological features.
- Current KO descriptions and structured KEGG evidence used in the final treatment-blind audit were snapshotted on 2026-10-08 and are preserved in the 7,931-KO review table.

## Utility-mapping limitation

The precise release, download date, and checksum of the original `map_ko_uniref90.txt.gz` file could not be recovered from the local analysis tree. No release identifier has been inferred. The reproducible analyses in this repository begin with the resulting KO relative-abundance matrices.

Database files themselves are not redistributed unless their licenses permit it.
