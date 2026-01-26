from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

def mol_to_svg(mol: Chem.Mol, add_h: bool = False, add_indices: bool = False, w: int = 300, h: int = 200) -> str:
    m = Chem.AddHs(mol) if add_h else mol
    
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m, canonOrient=True)

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.drawOptions().addAtomIndices = add_indices
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def run_rdkit(smiles: str, smarts: str) -> dict:
    """Run RDKit to get molecule and substructure match information."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Invalid SMILES input: '{smiles}'"}

    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return {"error": f"Invalid SMARTS input: '{smarts}'"}

    output = {}

    # Add canonical SMILES and InChI
    can = Chem.MolToSmiles(mol, canonical=True)
    inchi = Chem.MolToInchi(mol)
    can_from_inchi = Chem.MolToSmiles(Chem.MolFromInchi(inchi), canonical=True) if inchi else ""
    output["output_can"] = can
    output["output_inchi"] = inchi
    output["output_can_from_inchi"] = can_from_inchi

    # Add drawing SVG: plain, with hydrogens, and with indices
    output["svg"] = mol_to_svg(mol)
    output["svg_h"] = mol_to_svg(mol, add_h=True)
    output["svg_index"] = mol_to_svg(mol, add_indices=True)
    output["svg_h_index"] = mol_to_svg(mol, add_h=True, add_indices=True)

    # Find substructure matches
    matches = mol.GetSubstructMatches(patt)
    n = len(matches)
    if n == 0:
        output["n_matches"] = f"No matches found for SMARTS '{smarts}'."
    else:
        output["n_matches"] = f"Found {n} match{'es' if n > 1 else ''} for SMARTS '{smarts}'."
        msgs = [f" Matches at atom indices (0-based):"]
        msgs.append(" " + "\n ".join(str(match) for match in matches))
        output["smarts_matches"] = "\n".join(msgs)

    return output
