from openbabel import openbabel as ob

def mol_to_svg(mol: ob.OBMol, add_h: bool = False) -> str:
    m = ob.OBMol(mol)
    if add_h:
        m.AddHydrogens()

    op = ob.OBOp.FindType("gen2D")
    if op is not None:
        op.Do(m)

    conv = ob.OBConversion()
    if not conv.SetOutFormat("svg"):
        return ""
    return conv.WriteString(m) or ""

def _ob_write(mol: ob.OBMol, format: str) -> str:
    conv = ob.OBConversion()
    if not conv.SetOutFormat(format):
        return ""
    return conv.WriteString(mol) or ""

def _ob_read(infmt: str, s: str) -> ob.OBMol | None:
    conv = ob.OBConversion()
    if not conv.SetInFormat(infmt):
        return None
    m = ob.OBMol()
    return m if conv.ReadString(m, s) else None

def run_openbabel(smiles: str, smarts: str) -> dict:
    """Run OpenBabel to get molecule and substructure match information."""
    
    mol = _ob_read("smi", smiles)
    if mol is None:
        return {"error": f"Invalid SMILES input: '{smiles}'"}

    output = {}

    # Add canonical SMILES and InChI
    output["output_can"] = _ob_write(mol, "can")
    output["output_inchi"] = _ob_write(mol, "inchi")
    mol_from_inchi = _ob_read("inchi", output["output_inchi"]) if output["output_inchi"] else None
    output["output_can_from_inchi"] = _ob_write(mol_from_inchi, "can") if mol_from_inchi else ""

    # Add drawing SVG: plain, with hydrogens, and with indices
    output["svg"] = mol_to_svg(mol)
    output["svg_h"] = mol_to_svg(mol, add_h=True)

    # SMARTS matching
    patt = ob.OBSmartsPattern()
    if not patt.Init(smarts):
        return output

    # Find substructure matches
    matches = []
    if patt.Match(mol):
        match_list = patt.GetUMapList()
        for match in match_list:
            matches.append(tuple(match))
    n = len(matches)
    if n == 0:
        output["n_matches"] = f"No matches found for SMARTS '{smarts}'."
    else:
        output["n_matches"] = f"Found {n} match{'es' if n > 1 else ''} for SMARTS '{smarts}'."
        msgs = [f" Matches at atom indices (1-based):"]
        msgs.append(" " + "\n ".join(str(match) for match in matches))
        output["smarts_matches"] = "\n".join(msgs)

    return output
