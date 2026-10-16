#!/usr/bin/env python3
"""Validate the portable publication repository and its reproduced outputs."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
VALIDATION.mkdir(parents=True, exist_ok=True)
RESULTS = ROOT / "results" / "reproduced"

checks: list[dict] = []


def record(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def require(name: str, condition: bool, detail: str) -> None:
    record(name, bool(condition), detail)


def close(value: float, expected: float, tolerance: float = 5e-4) -> bool:
    return math.isfinite(value) and abs(value - expected) <= tolerance


def read_csv(relative: str, **kwargs) -> pd.DataFrame:
    path = ROOT / relative
    require(f"file exists: {relative}", path.is_file() and path.stat().st_size > 0,
            f"{path.stat().st_size if path.exists() else 0} bytes")
    return pd.read_csv(path, **kwargs)


# The package must be portable and self-contained.
symlinks = [str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_symlink()]
require("no symbolic links", not symlinks, "; ".join(symlinks) or "none")

forbidden = []
# Assemble machine-specific tokens so this validator does not flag its own source.
forbidden_tokens = ("/home/" + "felius", "wsl" + ".localhost")
for folder in (ROOT / "scripts", ROOT / "docs"):
    if not folder.exists():
        continue
    for path in folder.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".r", ".py", ".sh", ".md", ".txt", ".csv", ".tsv"}:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            if any(token in text for token in forbidden_tokens):
                forbidden.append(str(path.relative_to(ROOT)))
require("no machine-specific absolute paths", not forbidden, "; ".join(forbidden) or "none")

# Frozen inputs and curation invariants.
mapping = read_csv("data/functional_metagenomics/frozen_332_KO_functional_group_mapping.csv")
definitions = read_csv("data/functional_metagenomics/functional_group_definitions.csv")
universe = read_csv("data/functional_metagenomics/full_7931_KO_universe_review.csv")
require("332 unique curated KOs", len(mapping) == 332 and mapping.ko_id.nunique() == 332,
        f"rows={len(mapping)}, unique={mapping.ko_id.nunique()}")
require("13 functional groups", mapping.functional_group.nunique() == 13 and len(definitions) == 13,
        f"mapping={mapping.functional_group.nunique()}, definitions={len(definitions)}")
require("one primary group per KO", not mapping.ko_id.duplicated().any(), "no duplicate ko_id")
require("7,931-KO frozen universe", len(universe) == 7931 and universe.ko_id.nunique() == 7931,
        f"rows={len(universe)}, unique={universe.ko_id.nunique()}")
require("curated KOs belong to frozen universe", set(mapping.ko_id) <= set(universe.ko_id),
        f"missing={len(set(mapping.ko_id)-set(universe.ko_id))}")

# Final sample sizes from the reproduced analyses.
amp = read_csv("results/reproduced/amplicons/amplicon_global_statistics.csv")
expected_amp = {("field", "ITS"): 18, ("field", "16S"): 19,
                ("microcosm", "ITS"): 40, ("microcosm", "16S"): 36}
observed_amp = {(r.experiment, r.marker): int(r.N) for r in amp.itertuples()}
require("amplicon final sample sizes", observed_amp == expected_amp,
        f"observed={observed_amp}")

fun = read_csv("results/reproduced/functional_metagenomics/global_statistics_all_experiments.csv")
observed_fun = dict(zip(fun.experiment, fun.N.astype(int)))
require("functional metagenomics final sample sizes", observed_fun == {"field": 19, "microcosm": 41},
        f"observed={observed_fun}")
require("field analyzes 13 variable groups", int(fun.loc[fun.experiment.eq("field"), "analyzed_groups"].iloc[0]) == 13,
        str(fun.loc[fun.experiment.eq("field"), "analyzed_groups"].iloc[0]))
require("microcosm analyzes 12 variable groups", int(fun.loc[fun.experiment.eq("microcosm"), "analyzed_groups"].iloc[0]) == 12,
        str(fun.loc[fun.experiment.eq("microcosm"), "analyzed_groups"].iloc[0]))

# Numerical anchors reported in the final audit.
anchors = {
    ("field", "ITS"): (1.0993, 0.12784, 0.2953, 0.6573),
    ("field", "16S"): (1.18977, 0.12947, 0.1139, 0.1626),
    ("microcosm", "ITS"): (12.7476, 0.40795, 0.0001, 0.1282),
    ("microcosm", "16S"): (4.47695, 0.21342, 0.0001, 0.0230),
}
for key, target in anchors.items():
    row = amp[(amp.experiment == key[0]) & (amp.marker == key[1])].iloc[0]
    values = (row.F, row.R2, row.p, row.dispersion_p)
    passed = all(close(float(v), float(t), 7e-4) for v, t in zip(values, target))
    require(f"amplicon numerical anchor: {key[0]} {key[1]}", passed,
            f"F={row.F:.5g}, R2={row.R2:.5g}, p={row.p:.5g}, dispersion_p={row.dispersion_p:.5g}")

fun_targets = {
    "field": (2.39068, 0.23008, 0.0159, 0.0233),
    "microcosm": (4.1345, 0.178716, 0.0010, 0.3060),
}
for experiment, target in fun_targets.items():
    row = fun[fun.experiment.eq(experiment)].iloc[0]
    values = (row.F, row.R2, row.p, row.dispersion_p)
    passed = all(close(float(v), float(t), 8e-4) for v, t in zip(values, target))
    require(f"functional numerical anchor: {experiment}", passed,
            f"F={row.F:.5g}, R2={row.R2:.5g}, p={row.p:.5g}, dispersion_p={row.dispersion_p:.5g}")

measurements = read_csv("results/reproduced/functional_measurements/omnibus_tests.csv")
require("four functional-measurement endpoints", len(measurements) == 4,
        f"rows={len(measurements)}")
require("functional-measurement N values", set(measurements.N.astype(int)) == {43, 46},
        f"N={sorted(measurements.N.astype(int).unique())}")

table1 = read_csv("results/reproduced/chemistry/Table_1_descriptive_no_inference.csv")
require("three descriptive substrate rows", len(table1) == 3, f"rows={len(table1)}")

s1 = read_csv("tables/Table_S1_functional_group_definitions.csv")
s2 = read_csv("tables/Table_S2_complete_KO_mapping.csv")
s3 = read_csv("tables/Table_S3_functional_group_sensitivity.csv")
require("Table S1 has 13 groups", len(s1) == 13, f"rows={len(s1)}")
require("Table S2 has 332 KOs", len(s2) == 332 and s2.ko_id.nunique() == 332,
        f"rows={len(s2)}, unique={s2.ko_id.nunique()}")
require("Table S3 has 25 analyzed group/experiment rows", len(s3) == 25, f"rows={len(s3)}")

# Final figure assets are intentionally not regenerated by run_all.sh.
figure_names = ["Figure1.pptx"] + [f"Figure{i}.pdf" for i in range(2, 7)] + [f"FigS{i}.pdf" for i in range(1, 11)]
missing_figures = []
for name in figure_names:
    path = ROOT / "figures" / "final" / name
    if not path.is_file() or path.stat().st_size == 0:
        missing_figures.append(name)
require("final manuscript figure assets present", not missing_figures,
        ", ".join(missing_figures) or "Figures 1-6 and S1-S10 present")

# Basic PDF/PNG integrity signatures.
bad_signatures = []
for path in (ROOT / "figures" / "final").iterdir():
    if not path.is_file():
        continue
    head = path.read_bytes()[:8]
    if path.suffix.lower() == ".pdf" and not head.startswith(b"%PDF"):
        bad_signatures.append(path.name)
    if path.suffix.lower() == ".png" and head != b"\x89PNG\r\n\x1a\n":
        bad_signatures.append(path.name)
require("figure file signatures", not bad_signatures, ", ".join(bad_signatures) or "valid")

failed = [item for item in checks if item["status"] == "FAIL"]
report_lines = [
    "Publication repository validation",
    "=================================",
    "Repository: publication_repository (resolved relative to this script)",
    f"Checks: {len(checks)}; passed: {len(checks)-len(failed)}; failed: {len(failed)}",
    "",
]
report_lines.extend(f"[{item['status']}] {item['check']}: {item['detail']}" for item in checks)
(VALIDATION / "VALIDATION_REPORT.txt").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
(VALIDATION / "VALIDATION_REPORT.json").write_text(
    json.dumps({"passed": not failed, "checks": checks}, indent=2) + "\n", encoding="utf-8"
)

# Hash every release file except the manifest itself and ephemeral caches.
manifest_rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or path.name == "MANIFEST.sha256" or "__pycache__" in path.parts:
        continue
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_rows.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
(VALIDATION / "MANIFEST.sha256").write_text("\n".join(manifest_rows) + "\n", encoding="utf-8")

print("\n".join(report_lines[:5]))
if failed:
    for item in failed:
        print(f"FAIL: {item['check']}: {item['detail']}", file=sys.stderr)
    raise SystemExit(1)
print("All release checks passed.")
