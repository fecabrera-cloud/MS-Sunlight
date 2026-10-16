#!/usr/bin/env python3
"""Build publication-readable Tables S1-S3 and the complete KO audit workbook."""

from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.shared import Inches, Mm, Pt, RGBColor
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "reproduced" / "functional_metagenomics"
TABLES = ROOT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

DARK = "1F4E78"
LIGHT = "EEF5FA"
MID = "D9EAF7"
WHITE = "FFFFFF"


def make_s3(source: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in source.itertuples(index=False):
        primary = "Significant" if row.primary_significant else "Not significant"
        if row.primary_significant:
            if row.log_conclusion_stable and row.loo_conclusion_stable:
                interpretation = "Robust to log transformation and all single-sample omissions"
            elif row.log_conclusion_stable:
                interpretation = "Stable after log transformation; sensitive to single-sample omission"
            elif row.loo_conclusion_stable:
                interpretation = "Sensitive to log transformation; stable to single-sample omission"
            else:
                interpretation = "Sensitive to transformation and single-sample omission"
        elif row.log_significant or not row.loo_conclusion_stable:
            interpretation = (
                "Primary result nonsignificant; significance threshold crossed after "
                "some sensitivity analyses"
            )
        else:
            interpretation = "Stable nonsignificant result"
        rows.append(
            {
                "Experiment": row.experiment.title(),
                "Functional group": row.functional_group,
                "Primary BH-adjusted p": row.primary_q,
                "Log-transformed BH-adjusted p": row.log_q,
                "Leave-one-out minimum BH-adjusted p": row.loo_q_min,
                "Leave-one-out maximum BH-adjusted p": row.loo_q_max,
                "Omissions changing significance (n)": row.n_omissions_changing_significance,
                "Omissions changing significance (sample IDs)": row.omissions_changing_significance,
                "Primary result": primary,
                "Sensitivity interpretation": interpretation,
            }
        )
    return pd.DataFrame(rows)


def shade(cell, fill: str) -> None:
    props = cell._tc.get_or_add_tcPr()
    node = props.find(qn("w:shd"))
    if node is None:
        node = OxmlElement("w:shd")
        props.append(node)
    node.set(qn("w:fill"), fill)


def repeat_header(row) -> None:
    props = row._tr.get_or_add_trPr()
    node = OxmlElement("w:tblHeader")
    node.set(qn("w:val"), "true")
    props.append(node)


def no_split(row) -> None:
    props = row._tr.get_or_add_trPr()
    props.append(OxmlElement("w:cantSplit"))


def cell_text(cell, value, size=6.4, bold=False, color="000000", center=False) -> None:
    text = "—" if pd.isna(value) else str(value)
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.paragraph_format.space_after = Pt(0)
    paragraph.paragraph_format.line_spacing = 1
    if center:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run(text)
    run.font.name = "Arial"
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = RGBColor.from_string(color)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def hyperlink(cell, text: str, url: str) -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    relation = paragraph.part.relate_to(url, RT.HYPERLINK, is_external=True)
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), relation)
    run = OxmlElement("w:r")
    props = OxmlElement("w:rPr")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "0563C1")
    underline = OxmlElement("w:u")
    underline.set(qn("w:val"), "single")
    size = OxmlElement("w:sz")
    size.set(qn("w:val"), "13")
    props.extend([color, underline, size])
    run.append(props)
    label = OxmlElement("w:t")
    label.text = text
    run.append(label)
    link.append(run)
    paragraph._p.append(link)


def fixed_table(table, widths) -> None:
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    table._tbl.tblPr.append(layout)
    for row in table.rows:
        for index, width in enumerate(widths):
            if index < len(row.cells):
                row.cells[index].width = Inches(width)


def header(row) -> None:
    for cell in row.cells:
        shade(cell, DARK)
        for run in cell.paragraphs[0].runs:
            run.font.color.rgb = RGBColor(255, 255, 255)
            run.font.bold = True
    repeat_header(row)
    no_split(row)


def caption(document, text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.paragraph_format.keep_with_next = True
    run = paragraph.add_run(text)
    run.font.name = "Arial"
    run.font.size = Pt(9)
    run.font.bold = True


def build_word(s1: pd.DataFrame, s2: pd.DataFrame, s3: pd.DataFrame) -> None:
    document = Document()
    section = document.sections[0]
    section.orientation = WD_ORIENT.LANDSCAPE
    section.page_width, section.page_height = Mm(297), Mm(210)
    for side in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(section, side, Inches(0.45))
    document.styles["Normal"].font.name = "Arial"
    document.styles["Normal"].font.size = Pt(8)
    title = document.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title_run = title.add_run("Supplementary Tables S1–S3")
    title_run.font.name = "Arial"
    title_run.font.size = Pt(14)
    title_run.font.bold = True

    caption(
        document,
        "Table S1. Definitions and membership of the 13 targeted functional groups used "
        "in the functional metagenomic analysis. The shared KO mapping was applied to "
        "both experiments; lignin-associated oxidation was not detected in the final "
        "microcosm dataset, leaving 12 analyzable groups.",
    )
    keys = [
        "display_order", "functional_group", "hypothesis_role", "process_represented",
        "inclusion_criterion", "important_exclusions", "n_kos_shared_mapping",
        "n_kos_detected_field", "n_kos_detected_microcosm",
    ]
    labels = [
        "Order", "Functional group", "Hypothesis role", "Process represented",
        "Inclusion criterion", "Important exclusions", "KOs in mapping",
        "Detected field", "Detected microcosm",
    ]
    table = document.add_table(rows=1, cols=len(keys))
    for index, label in enumerate(labels):
        cell_text(table.rows[0].cells[index], label, bold=True, color=WHITE, center=True)
    header(table.rows[0])
    for ridx, row in s1.sort_values("display_order").iterrows():
        cells = table.add_row().cells
        for index, key in enumerate(keys):
            cell_text(cells[index], row[key], center=index in {0, 6, 7, 8})
            if ridx % 2:
                shade(cells[index], LIGHT)
        no_split(table.rows[-1])
    fixed_table(table, [0.4, 1.4, 1.0, 1.55, 1.75, 1.55, 0.7, 0.7, 0.8])

    document.add_page_break()
    caption(
        document,
        "Table S2. Manually curated KEGG Orthology assignments to targeted functional "
        "groups. Each KO was assigned to one primary functional group to prevent double "
        "counting. This Word table is an abridged presentation; the accompanying "
        "Table_S2_complete_KO_mapping.xlsx contains the complete audit record.",
    )
    labels = ["KO", "KEGG annotation", "Functional group", "Confidence", "Curation justification", "Detection"]
    table = document.add_table(rows=1, cols=6)
    for index, label in enumerate(labels):
        cell_text(table.rows[0].cells[index], label, bold=True, color=WHITE, center=True)
    header(table.rows[0])
    order = dict(zip(s1.functional_group, s1.display_order))
    s2 = s2.assign(_order=s2.functional_group.map(order)).sort_values(["_order", "ko_id"])
    last = None
    for _, row in s2.iterrows():
        if row.functional_group != last:
            group_row = table.add_row()
            merged = group_row.cells[0].merge(group_row.cells[-1])
            count = int((s2.functional_group == row.functional_group).sum())
            cell_text(merged, f"{row.functional_group} (n = {count} KOs)", size=7, bold=True)
            shade(merged, MID)
            no_split(group_row)
            last = row.functional_group
        cells = table.add_row().cells
        hyperlink(cells[0], row.ko_id, row.kegg_entry_url)
        detection = (
            "Field and microcosm" if row.detected_field and row.detected_microcosm
            else "Field only" if row.detected_field
            else "Microcosm only" if row.detected_microcosm
            else "Not detected"
        )
        for index, value in enumerate(
            [row.kegg_name, row.functional_group, row.assignment_confidence, row.curation_reason, detection],
            start=1,
        ):
            cell_text(cells[index], value, center=index in {3, 5})
        no_split(table.rows[-1])
    fixed_table(table, [0.7, 2.2, 1.45, 0.75, 4.1, 1.1])

    document.add_page_break()
    caption(
        document,
        "Table S3. Sensitivity of treatment effects on targeted functional-group relative "
        "abundances. Primary BH-adjusted p-values are shown with natural-log-transformation "
        "and leave-one-sample-out results. These diagnostics did not trigger sample exclusion "
        "or replace the primary full-data analysis.",
    )
    table = document.add_table(rows=1, cols=len(s3.columns))
    short = ["Experiment", "Functional group", "Primary q", "Log-scale q", "LOO q min", "LOO q max", "Changed n", "Sample IDs", "Primary result", "Interpretation"]
    for index, label in enumerate(short):
        cell_text(table.rows[0].cells[index], label, bold=True, color=WHITE, center=True)
    header(table.rows[0])
    for ridx, row in s3.iterrows():
        cells = table.add_row().cells
        for index, value in enumerate(row):
            if isinstance(value, float) and index in {2, 3, 4, 5}:
                value = f"{value:.4f}"
            cell_text(cells[index], value, size=6.1, center=index in {0, 2, 3, 4, 5, 6})
            if ridx % 2:
                shade(cells[index], LIGHT)
        no_split(table.rows[-1])
    fixed_table(table, [0.7, 1.4, 0.65, 0.68, 0.68, 0.68, 0.65, 1.25, 1.0, 2.6])
    document.save(TABLES / "Supplementary_Tables_S1-S3.docx")


def build_excel(s2: pd.DataFrame) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Complete KO mapping"
    columns = list(s2.columns)
    for index, name in enumerate(columns, start=1):
        cell = sheet.cell(1, index, name)
        cell.fill = PatternFill("solid", fgColor=DARK)
        cell.font = Font(name="Arial", bold=True, color=WHITE)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    for ridx, row in enumerate(s2.itertuples(index=False, name=None), start=2):
        for cidx, value in enumerate(row, start=1):
            sheet.cell(ridx, cidx, None if pd.isna(value) else value)
        ko = sheet.cell(ridx, columns.index("ko_id") + 1)
        ko.hyperlink = s2.iloc[ridx - 2].kegg_entry_url
        ko.font = Font(name="Arial", color="0563C1", underline="single")
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = f"A1:{get_column_letter(len(columns))}{len(s2)+1}"
    for index, name in enumerate(columns, start=1):
        sheet.column_dimensions[get_column_letter(index)].width = min(
            70, max(12, 2 + max(len(name), *(len(str(v)) for v in s2[name].dropna().head(200))))
        )
    documentation = workbook.create_sheet("Documentation")
    notes = [
        ("Rows", len(s2)),
        ("Unique KO identifiers", s2.ko_id.nunique()),
        ("Functional groups", s2.functional_group.nunique()),
        ("Assignment rule", "Each KO is assigned to one primary group to prevent double counting."),
        ("Scope", "Targeted hypothesis-driven mapping, not an exhaustive metabolic profile."),
    ]
    for index, (field, value) in enumerate(notes, start=1):
        documentation.cell(index, 1, field).font = Font(bold=True, color=DARK)
        documentation.cell(index, 2, value)
    documentation.column_dimensions["A"].width = 28
    documentation.column_dimensions["B"].width = 100
    workbook.save(TABLES / "Table_S2_complete_KO_mapping.xlsx")


def main() -> None:
    s1 = pd.read_csv(RESULTS / "functional_group_definitions.csv")
    s2 = pd.read_csv(RESULTS / "frozen_332_KO_functional_group_mapping.csv")
    sensitivity = pd.read_csv(RESULTS / "functional_group_sensitivity_summary.csv")
    s3 = make_s3(sensitivity)
    s1.to_csv(TABLES / "Table_S1_functional_group_definitions.csv", index=False)
    s2.to_csv(TABLES / "Table_S2_complete_KO_mapping.csv", index=False)
    s3.to_csv(TABLES / "Table_S3_functional_group_sensitivity.csv", index=False)
    s3.to_excel(TABLES / "Table_S3_functional_group_sensitivity.xlsx", index=False)
    build_word(s1, s2, s3)
    build_excel(s2)


if __name__ == "__main__":
    main()
