import os
import pubchempy
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Draw
from flask import Flask, render_template, request, send_from_directory, jsonify
import requests
from bs4 import BeautifulSoup
from io import BytesIO
from Bio import Entrez
import logging
from openpyxl import Workbook
from openpyxl.drawing.image import Image as ExcelImage
from datetime import datetime
from time import sleep
from random import randint

# Set email for Entrez
Entrez.email = "your_email@example.com"

# Flask app setup
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'static/cadd_images'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Set up logging
logging.basicConfig(level=logging.INFO)

# --------------------------------
# Utility Functions
# --------------------------------

def fetch_pubchem_data(name):
    """Fetch compound data from PubChem."""
    try:
        cids = pubchempy.get_cids(name, namespace="name")
        if not cids:
            return None
        compound = pubchempy.get_compounds(cids[0], namespace="cid")
        return compound[0] if compound else None
    except Exception as e:
        logging.error(f"Error fetching PubChem data for {name}: {e}")
        return None


def generate_cadd_image(smiles, filename):
    """Generate CADD image from SMILES and save it."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        img = Draw.MolToImage(mol)
        img.save(filename)
        return filename
    except Exception as e:
        logging.error(f"Error generating CADD image: {e}")
        return None


def determine_interaction(binding_affinity1, binding_affinity2, synergy_binding_affinity):
    """Determine interaction type based on binding affinities."""
    if synergy_binding_affinity < (binding_affinity1 + binding_affinity2):
        return "Agonistic", "Improved binding affinity."
    elif synergy_binding_affinity > (binding_affinity1 + binding_affinity2):
        return "Antagonistic", "Reduced binding affinity."
    return "Neutral", "No significant change."


@app.route('/')
def home():
    return render_template("index.html")

def save_results_with_images(results, file_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "Synergy Results"

    # Define headers
    headers = [
        "Ingredient 1", "Ingredient 2", "Binding Affinity",
        "Interaction", "Remarks", "Individual Effects", 
        "Synergy Effects", "PubMed Data (1)", "PubMed Data (2)", 
        "Synergy PubMed Data", "CADD Analytics Image (1)", 
        "CADD Analytics Image (2)", "Synergy Image"
    ]
    ws.append(headers)

    for row_idx, result in enumerate(results, start=2):
        ws.append([
            result["Ingredient 1"], result["Ingredient 2"], result["Binding Affinity"],
            result["Interaction"], result["Remarks"], result["Individual Effects"],
            result["Synergy Effects"], ", ".join(
                [f"{pub['PubMed ID']}: {pub['Title']}" for pub in result["PubMed Data 1"]]
            ), ", ".join(
                [f"{pub['PubMed ID']}: {pub['Title']}" for pub in result["PubMed Data 2"]]
            ), ", ".join(
                [f"{pub['PubMed ID']}: {pub['Title']}" for pub in result["Synergy PubMed Data"]]
            ), None, None, None  # Placeholder for images
        ])

        # Embed images
        if result["CADD Analytics Image 1"]:
            img1 = ExcelImage(result["CADD Analytics Image 1"])
            img1.anchor = f"L{row_idx}"
            ws.add_image(img1)

        if result["CADD Analytics Image 2"]:
            img2 = ExcelImage(result["CADD Analytics Image 2"])
            img2.anchor = f"M{row_idx}"
            ws.add_image(img2)

        if result["Synergy Image"]:
            img3 = ExcelImage(result["Synergy Image"])
            img3.anchor = f"N{row_idx}"
            ws.add_image(img3)

    wb.save(file_path)

def fetch_pubmed_ids_titles(query):
    """Fetch PubMed IDs and Titles based on a query."""
    try:
        search_results = Entrez.esearch(db="pubmed", term=query, retmax=5)
        record = Entrez.read(search_results)
        ids = record.get("IdList", [])
        articles = []

        for pubmed_id in ids:
            fetch_result = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="abstract", retmode="text")
            article_data = fetch_result.read()
            articles.append({"PubMed ID": pubmed_id, "Title": article_data.split("\n")[0]})
        
        return articles
    except Exception as e:
        logging.error(f"Error fetching PubMed data for query {query}: {e}")
        return []



@app.route('/cadd', methods=['GET', 'POST'])
def cadd_synergy():
    if request.method == 'POST':
        ingredient1 = request.form['ingredient1']
        ingredient2 = request.form['ingredient2']
        results = []

        # Fetch compound data
        compound1 = fetch_pubchem_data(ingredient1)
        compound2 = fetch_pubchem_data(ingredient2)

        if not compound1 or not compound2:
            return render_template('cadd.html', error="Failed to fetch data for one or both ingredients.")

        # Generate CADD images
        cadd1_img = generate_cadd_image(compound1.canonical_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient1}_CADD.png"))
        cadd2_img = generate_cadd_image(compound2.canonical_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient2}_CADD.png"))
        synergy_smiles = f"{compound1.canonical_smiles}.{compound2.canonical_smiles}"
        synergy_cadd_img = generate_cadd_image(synergy_smiles, os.path.join(app.config['UPLOAD_FOLDER'], f"{ingredient1}_{ingredient2}_Synergy.png"))

        # Dummy binding affinities
        binding_affinity1 = 10.0
        binding_affinity2 = 10.0
        synergy_binding_affinity = (binding_affinity1 + binding_affinity2) * 0.9

        # Determine interaction
        interaction, remarks = determine_interaction(binding_affinity1, binding_affinity2, synergy_binding_affinity)

        # Fetch PubMed data
        pubmed_data1 = fetch_pubmed_ids_titles(ingredient1)
        pubmed_data2 = fetch_pubmed_ids_titles(ingredient2)
        synergy_pubmed_data = fetch_pubmed_ids_titles(f"{ingredient1} {ingredient2}")

        # Prepare results
        results.append({
            "Ingredient 1": ingredient1,
            "Ingredient 2": ingredient2,
            "Binding Affinity": f"{binding_affinity1}, {binding_affinity2}",
            "Interaction": interaction,
            "Remarks": remarks,
            "CADD Analytics Image 1": cadd1_img,
            "CADD Analytics Image 2": cadd2_img,
            "Synergy Image": synergy_cadd_img,
            "PubMed Data 1": pubmed_data1,
            "PubMed Data 2": pubmed_data2,
            "Synergy PubMed Data": synergy_pubmed_data,
            "Individual Effects": "Improves metabolism and reduces inflammation.",  # Example effect
            "Synergy Effects": "Enhances anti-inflammatory effects and reduces toxicity.",  # Example synergy effect
        })

        # Export results to Excel with images
        excel_file = os.path.join(app.config['UPLOAD_FOLDER'], "synergy_results.xlsx")
        save_results_with_images(results, excel_file)
        return render_template('cadd.html', results=results, download_link="synergy_results.xlsx")

    return render_template('cadd.html')


@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename, as_attachment=True)

# Default query for PubMed
DEFAULT_QUERY = "((((((((Vitamin) OR (mineral)) OR (nutraceutical)) OR (phytochemical)) OR (plants)) OR (herb)) OR (botanical)) OR (health supplement)) OR (medicinal plant) AND ((classicalarticle[Filter] OR clinicalstudy[Filter] OR clinicaltrial[Filter] OR governmentpublication[Filter] OR meta-analysis[Filter] OR multicenterstudy[Filter] OR observationalstudy[Filter] OR overall[Filter] OR practiceguideline[Filter] OR randomizedcontrolledtrial[Filter] OR researchsupportnonusgovt[Filter] OR researchsupportusgovtnonphs[Filter] OR researchsupportusgovtphs[Filter] OR researchsupportusgovernment[Filter] OR systematicreview[Filter]) AND (humans[Filter]) AND (english[Filter]))"

# Function to fetch metadata from DOI URL with retry logic
def fetch_doi_metadata(doi, headers):
    if not doi:
        return ["N/A"] * 7

    doi_url = f"https://doi.org/{doi}"
    retries = 3
    for attempt in range(retries):
        try:
            response = requests.get(doi_url, headers=headers, timeout=10)
            if response.status_code != 200:
                continue
            soup = BeautifulSoup(response.text, "html.parser")

            fields = ["citation_journal_title", "citation_volume", "citation_issue", 
                      "citation_firstpage", "citation_publisher", "citation_date", "citation_issn"]
            metadata = [
                (soup.find("meta", {"name": field})["content"] if soup.find("meta", {"name": field}) else "N/A")
                for field in fields
            ]
            return metadata
        except requests.RequestException as e:
            print(f"Error fetching DOI metadata (attempt {attempt + 1}): {e}")
            sleep(randint(1, 3))

    return ["N/A"] * 7

# Function to scrape PubMed articles
def scrape_pubmed_articles(query, start_date, end_date):
    base_url = "https://pubmed.ncbi.nlm.nih.gov/"
    headers = {"User-Agent": "Mozilla/5.0"}
    query = f"({query}) AND ({start_date.strftime('%Y/%m/%d')}:{end_date.strftime('%Y/%m/%d')}[pdat])"
    response = requests.get(base_url, params={"term": query, "size": 200}, headers=headers)
    soup = BeautifulSoup(response.text, "html.parser")
    articles = []

    for article in soup.find_all("article", class_="full-docsum"):
        try:
            title_elem = article.find("a", class_="docsum-title")
            title = title_elem.text.strip() if title_elem else "N/A"
            url = f"{base_url}{title_elem['href']}" if title_elem else "N/A"

            authors_elem = article.find("span", class_="docsum-authors")
            authors = authors_elem.text.strip() if authors_elem else "N/A"

            published_date_elem = article.find("span", class_="docsum-journal-citation")
            published_date = published_date_elem.text.strip().split(";")[0] if published_date_elem else "N/A"

            journal_elem = article.find("span", class_="docsum-journal-name")
            journal = journal_elem.text.strip() if journal_elem else "N/A"

            pmid_elem = article.find("span", class_="docsum-pmid")
            pmid = pmid_elem.text.strip() if pmid_elem else "N/A"

            doi_elem = article.find("span", class_="docsum-doi")
            doi = doi_elem.text.strip() if doi_elem else "N/A"

            abstract_elem = article.find("div", class_="docsum-abstract")
            abstract = abstract_elem.text.strip() if abstract_elem else "N/A"

            if doi != "N/A":
                journal, volume, issue, pages, publisher, date_published, issn = fetch_doi_metadata(doi, headers)
            else:
                volume, issue, pages, publisher, date_published, issn = ["N/A"] * 6

            publication_year = published_date.split()[-1] if published_date != "N/A" else "N/A"

            articles.append({
                "Item type": "Journal Article",
                "Authors": authors,
                "Title": title,
                "Journal": journal,
                "Full journal": f"{journal}, {volume}, {issue}, {pages}",
                "Published Date": published_date,
                "Publication year": publication_year,
                "Volume": volume,
                "Issue": issue,
                "Pages": pages,
                "Publisher": publisher,
                "Date published": date_published,
                "ISSN": issn,
                "PMID": pmid,
                "DOI": doi,
                "URLs": url,
                "Abstract": abstract
            })

        except Exception as e:
            logging.error(f"Error processing article: {e}")
            continue

    return articles

@app.route('/pubmed', methods=['GET', 'POST'])
def pubmed_scraper():
    if request.method == 'POST':
        try:
            start_date = datetime.strptime(request.form.get("start_date"), "%Y-%m-%d")
            end_date = datetime.strptime(request.form.get("end_date"), "%Y-%m-%d")
            use_default_query = request.form.get("use_default_query")

            query = DEFAULT_QUERY if use_default_query else request.form.get("query", "")
            articles = scrape_pubmed_articles(query, start_date, end_date)

            if articles:
                excel_file = os.path.join(app.config['UPLOAD_FOLDER'], "articles.xlsx")
                df = pd.DataFrame(articles)
                df.to_excel(excel_file, index=False)
                return render_template(
                    "pubmed.html",
                    articles=articles,
                    download_link="articles.xlsx",
                    query=query,
                    use_default_query=use_default_query,
                    DEFAULT_QUERY=DEFAULT_QUERY
                )

            return render_template(
                "pubmed.html",
                error="No articles found.",
                query=query,
                use_default_query=use_default_query,
                DEFAULT_QUERY=DEFAULT_QUERY
            )
        except Exception as e:
            logging.error(f"Error in pubmed_scraper: {e}")
            return render_template(
                "pubmed.html",
                error="An error occurred while processing your request.",
                query="",
                use_default_query=False,
                DEFAULT_QUERY=DEFAULT_QUERY
            )
    return render_template("pubmed.html", query="", use_default_query=False, DEFAULT_QUERY=DEFAULT_QUERY)









if __name__ == "__main__":
    app.run(debug=True)
