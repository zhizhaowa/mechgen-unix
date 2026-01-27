import argparse
from flask import Flask, render_template, request

from chem_rdkit import run_rdkit
from chem_ob import run_openbabel

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    smiles = ""
    smarts = ""
    smirks = ""
    result = None

    if request.method == 'POST':
        smiles = request.form.get('smiles', '').strip()
        smarts = request.form.get('smarts', '').strip()
        smirks = request.form.get('smirks', '').strip()
        rd = run_rdkit(smiles, smarts, smirks)
        ob = run_openbabel(smiles, smarts)
        result = {"rdkit": rd, "openbabel": ob}

    return render_template('index.html', smiles=smiles, smarts=smarts, smirks=smirks, result=result)

@app.route("/rdk", methods=["GET", "POST"])
def rdk():
    smiles = smarts = smirks = ""
    r = None

    if request.method == "POST":
        smiles = (request.form.get("smiles") or "").strip()
        smarts = (request.form.get("smarts") or "").strip()
        smirks = (request.form.get("smirks") or "").strip()
        r = run_rdkit(smiles, smarts, smirks)

    return render_template(
        "rdk.html",
        smiles=smiles, smarts=smarts, smirks=smirks, r=r
    )

if __name__ == "__main__":
    #app.run(host="0.0.0.0", port=5005)
    # Localhost
    app.run(host="127.0.0.1", port=5005)
