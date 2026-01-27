# !/usr/bin/env python3

"""
Generates SVG images from SMILES using RDKit.
"""

from flask import Flask, request, Response
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from openbabel import openbabel as ob

app = Flask(__name__)

rdDepictor.SetPreferCoordGen(True)


@app.get("/svg")
def svg():
    """Generates an SVG image from a SMILES string."""

    smiles = request.args.get("smiles", "").strip()
    with_h = request.args.get("w_h", "0") == "1"
    use_ob = True

    if not smiles:
        return Response("Missing 'smiles' parameter.", status=400)

    if use_ob:  # OpenBabel version
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "svg")

        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return Response(f"Invalid SMILES input: '{smiles}'", status=400)

        if with_h:
            mol.AddHydrogens()

        svg_text = conv.WriteString(mol)

    else:  # RDKit version
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return Response(f"Invalid SMILES input: '{smiles}'", status=400)

        if with_h:
            mol = Chem.AddHs(mol)

        rdDepictor.Compute2DCoords(mol, canonOrient=True)
        drawer = rdMolDraw2D.MolDraw2DSVG(360, 260)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_text = drawer.GetDrawingText()

    if not svg_text:
        return Response("Failed to generate SVG.", status=500)
    return Response(svg_text, mimetype="image/svg+xml")


if __name__ == "__main__":
    # app.run(host="0.0.0.0", port=5000)
    app.run(host="127.0.0.1", port=5000)
