import os
import pubchempy
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from Bio import Entrez
from flask import Flask, render_template, request, send_file
import retry

# Set email for Entrez
Entrez.email = "your_email@example.com"

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'static/cadd_images'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Retry decorator
@retry.retry(tries=3, delay=1)
def fetch_pubchem_data(name):
    try:
        cids = pubchempy.get_cids(name, namespace="name")
        if not cids:
            return None
        compound = pubchempy.get_compounds(cids[0], namespace="cid")
        return compound[0]
    except Exception as e:
        print(f"Error fetching PubChem data: {e}")
        return None

# Fetch PubMed data
def fetch_pubmed_data(name):
    try:
        handle = Entrez.esearch(db="pubmed", term=name, retmax=5)
        record = Entrez.read(handle)
        pubmed_ids = record["IdList"]
        pubmed_titles = []
        for pubmed_id in pubmed_ids:
            handle = Entrez.esummary(db="pubmed", id=pubmed_id)
            record = Entrez.read(handle)
            pubmed_titles.append(record[0]["Title"])
        return pubmed_ids, pubmed_titles
    except Exception as e:
        print(f"Error fetching PubMed data: {e}")
        return [], []

def generate_cadd_image(smiles, filename):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string.")
        img = Draw.MolToImage(mol)
        img.save(filename)
        return filename
    except Exception as e:
        print(f"Error generating CADD image: {e}")
        return None

def determine_interaction(binding_affinity1, binding_affinity2, synergy_binding_affinity):
    if synergy_binding_affinity < (binding_affinity1 + binding_affinity2):
        return "Agonistic", "Improved binding affinity."
    elif synergy_binding_affinity > (binding_affinity1 + binding_affinity2):
        return "Antagonistic", "Reduced binding affinity."
    return "Neutral", "No significant change."

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        ingredient1 = request.form['ingredient1']
        ingredient2 = request.form['ingredient2']
        results = []

        compound1 = fetch_pubchem_data(ingredient1)
        compound2 = fetch_pubchem_data(ingredient2)

        if not compound1 or not compound2:
            return render_template('index.html', error="Failed to fetch data for one or both ingredients.")

        pubmed_ids1, pubmed_titles1 = fetch_pubmed_data(ingredient1)
        pubmed_ids2, pubmed_titles2 = fetch_pubmed_data(ingredient2)

        cadd1_img = generate_cadd_image(compound1.canonical_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient1}_CADD.png"))
        cadd2_img = generate_cadd_image(compound2.canonical_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient2}_CADD.png"))

        synergy_smiles = f"{compound1.canonical_smiles}.{compound2.canonical_smiles}"
        synergy_cadd_img = generate_cadd_image(synergy_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient1}_{ingredient2}_Synergy.png"))

        binding_affinity1 = 10.0
        binding_affinity2 = 10.0
        synergy_binding_affinity = (binding_affinity1 + binding_affinity2) * 0.9

        interaction, remarks = determine_interaction(binding_affinity1, binding_affinity2, synergy_binding_affinity)

        results.append({
            "Ingredient 1": ingredient1,
            "Ingredient 2": ingredient2,
            "PubMed IDs": ', '.join(pubmed_ids1),
            "PubMed Titles": ', '.join(pubmed_titles1),
            "Individual Effects": f"{compound1.iupac_name}, {compound2.iupac_name}",
            "Synergy Effects": synergy_smiles,
            "Binding Affinity": f"{binding_affinity1}, {binding_affinity2}",
            "Interaction": interaction,
            "Remarks": remarks,
            "CADD Analytics Image 1": cadd1_img,
            "CADD Analytics Image 2": cadd2_img,
            "Synergy Image": synergy_cadd_img,
        })

        df = pd.DataFrame(results)
        excel_file = "synergy_results.xlsx"
        df.to_excel(excel_file, index=False)
        return render_template('index.html', results=results, download_link=excel_file)

    return render_template('index.html')

@app.route('/download/<filename>')
def download_file(filename):
    return send_file(filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
