"""Assemble the numbered supplementary figures from the final PDF assets."""
from pathlib import Path
import json
import pypdfium2 as pdfium
from docx import Document
from docx.enum.section import WD_SECTION_START
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
captions = json.loads((ROOT / "docs/supplementary_captions.json").read_text())
out = ROOT / "results/reproduced/supplementary_figures"
out.mkdir(parents=True, exist_ok=True)
document = Document()
document.core_properties.title = "Supplementary Figures S1–S11"
document.core_properties.author = ""
document.core_properties.last_modified_by = ""
document.core_properties.comments = ""
document.styles["Normal"].font.name = "Arial"
document.styles["Normal"].font.size = Pt(10)
for i in range(1, 12):
    name = f"FigS{i}"
    pdf = pdfium.PdfDocument(ROOT / "figures/final" / f"{name}.pdf")
    assert len(pdf) == 1, name
    page = pdf[0]
    width, height = page.get_size()
    section = document.sections[0] if i == 1 else document.add_section(WD_SECTION_START.NEW_PAGE)
    section.page_width = Inches(max(8.3, width / 72 + 1))
    section.page_height = Inches(max(11.7, height / 72 + 2.2))
    for side in ["top_margin", "bottom_margin", "left_margin", "right_margin"]:
        setattr(section, side, Inches(.5))
    document.add_paragraph(f"Figure S{i}.", style="Heading 1")
    image = out / f"{name}.png"
    page.render(scale=2).to_pil().save(image)
    document.add_picture(str(image), width=Inches(width / 72), height=Inches(height / 72))
    document.add_paragraph(captions[name])
    pdf.close()
document.save(ROOT / "figures/final/Supplementary-Figures.docx")
print("Assembled Supplementary-Figures.docx from the 11 final PDFs.")
