# !/usr/bin/env python3

"""
Generates SVG images from SMILES using RDKit.
"""

from flask import Flask, request, Response
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

app = Flask(__name__)

@app.get("/svg")
def svg():
    """Generates an SVG image from a SMILES string."""

    smiles = request.args.get("smiles", "")  # Get SMILES
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to molecule
    if mol is None:
        return Response(f"Invalid SMILES input: '{smiles}'", status=400)

    rdDepictor.Compute2DCoords(mol)  # Compute 2D coordinates
    svg_txt = Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(360, 260),
        useSVG=True
    )
    return Response(svg_txt, mimetype="image/svg+xml")

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5005)

