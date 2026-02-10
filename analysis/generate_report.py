from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
import os

doc = SimpleDocTemplate("outputs/analysis/MD_report.pdf", pagesize=letter)
styles = getSampleStyleSheet()
elements = []

elements.append(Paragraph("Molecular Dynamics Report", styles["Title"]))
elements.append(Spacer(1, 20))

figures = [
    ("Protein RMSD", "outputs/analysis/rmsd.png"),
    ("Protein RMSF", "outputs/analysis/rmsf.png"),
    ("Radius of Gyration", "outputs/analysis/rg.png"),
    ("Free Energy Landscape", "outputs/analysis/fel.png"),
    ("Protein–Ligand Distance", "outputs/analysis/ligand_distance.png"),
]

for title, path in figures:
    if os.path.exists(path):
        elements.append(Paragraph(title, styles["Heading2"]))
        elements.append(Image(path, width=400, height=300))
        elements.append(Spacer(1, 20))

doc.build(elements)

print("PDF report generated: outputs/analysis/MD_report.pdf")

