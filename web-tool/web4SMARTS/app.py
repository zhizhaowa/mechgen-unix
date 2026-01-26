import argparse
from flask import Flask, render_template, request

from chem_rdkit import run_rdkit
from chem_ob import run_openbabel


def create_app() -> Flask:
    app = Flask(__name__)

    @app.route('/', methods=['GET', 'POST'])
    def index():
        smiles = ""
        smarts = ""
        result = None

        if request.method == 'POST':
            smiles = request.form.get('smiles', '').strip()
            smarts = request.form.get('smarts', '').strip()
            rd = run_rdkit(smiles, smarts)
            ob = run_openbabel(smiles, smarts)
            result = {"rdkit": rd, "openbabel": ob}

        return render_template('index.html', smiles=smiles, smarts=smarts, result=result)
    return app


def main():
    p = argparse.ArgumentParser()
    args = p.parse_args()

    app = create_app()
    #app.run(host="0.0.0.0", port=5005)

    # Localhost
    app.run(host="127.0.0.1", port=5005)


if __name__ == "__main__":
    main()

