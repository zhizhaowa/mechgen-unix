from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor, AllChem, inchi

def normalized(m):
    m2 = Chem.Mol(m)
    Chem.SanitizeMol(m2)
    return m2

def mol_to_svg(mol: Chem.Mol, add_h: bool = False, add_indices: bool = False, w: int = 300, h: int = 200) -> str:
    m = Chem.AddHs(mol) if add_h else mol
    
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m, canonOrient=True)

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.drawOptions().addAtomIndices = add_indices
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def run_rdkit(smiles: str, smarts: str, smirks: str) -> dict:
    """Run RDKit to get molecule and substructure match information."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Invalid SMILES input: '{smiles}'"}

    output = {}

     # Add canonical SMILES and InChI
    can = Chem.MolToSmiles(mol, canonical=True)
    #inc = Chem.MolToInchi(mol)
    inc = inchi.MolToInchi(mol)
    can_from_inchi = Chem.MolToSmiles(Chem.MolFromInchi(inc), canonical=True) if inc else ""
    output["output_can"] = can
    output["output_inchi"] = inc
    output["output_can_from_inchi"] = can_from_inchi

    # Add drawing SVG: plain, with hydrogens, and with indices
    output["svg"] = mol_to_svg(mol)
    output["svg_h"] = mol_to_svg(mol, add_h=True)
    output["svg_index"] = mol_to_svg(mol, add_indices=True)
    output["svg_h_index"] = mol_to_svg(mol, add_h=True, add_indices=True)
 
    # Process SMARTS if provided
    if smarts:
        patt = Chem.MolFromSmarts(smarts)
        if not patt:
            output["error"] = f"Invalid SMARTS input: '{smarts}'"
            return output

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

    # Process SMIRKS if provided
    if smirks:
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception as e:
            output["error"] = f"Invalid SMIRKS input: '{smirks}'"
            return output
        if not rxn:
            output["error"] = f"Invalid SMIRKS input: '{smirks}'"
            return output
        # Only allow single reactant reactions by now
        n_reactants = rxn.GetNumReactantTemplates()
        if n_reactants != 1:
            output["error"] = f"Only single reactant SMIRKS are supported. Provided SMIRKS has {n_reactants} reactants."
            return output

        ps = rxn.RunReactants((mol,))
        n_pd = sum(len(p) for p in ps)

        if n_pd == 0:
            output["n_products"] = f"No products generated."
        else:
            output["n_products"] = f"Generated {n_pd} product{'s' if n_pd > 1 else ''} for SMIRKS:\n  '{smirks}'."
            prod_smiles = []
            prod_svgs = []
            for pset in ps:
                for p in pset:
                    p = normalized(p)
                    psmi = Chem.MolToSmiles(p, canonical=True)
                    prod_smiles.append(psmi)
                    psvg = mol_to_svg(p, add_h=True)
                    prod_svgs.append(psvg)
            output["product_smiles"] = "\n".join(prod_smiles)
            output["product_svgs"] = prod_svgs

    return output
