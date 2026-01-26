# !/usr/bin/env python3

"""
Generates SVG images from SMILES using RDKit.
"""

from flask import Flask, request, Response
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

app = Flask(__name__)

rdDepictor.SetPreferCoordGen(True)

@app.get("/svg")
def svg():
    """Generates an SVG image from a SMILES string."""

    smiles = request.args.get("smiles", "")  # Get SMILES
    with_h = request.args.get("w_h", "0") == "1"

    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to molecule
    if mol is None:
        return Response(f"Invalid SMILES input: '{smiles}'", status=400)

    if with_h:
        mol = Chem.AddHs(mol)

    # rdDepictor.Compute2DCoords(mol)  # Compute 2D coordinates
    # AllChem.Compute2DCoords(mol)

    # Single molecule
    # svg_txt = Draw.MolToSVG(mol, 360, 260)
    # Multiple molecules
    # svg_txt = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(360, 260), useSVG=True)

    rdDepictor.Compute2DCoords(mol, canonOrient=True)
    drawer = rdMolDraw2D.MolDraw2DSVG(360, 260)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    #return Response(svg_txt, mimetype="image/svg+xml")
    return Response(drawer.GetDrawingText(), mimetype="image/svg+xml")

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
    # app.run(host="127.0.0.1", post=5000)
