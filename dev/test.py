import requests
import time
from bs4 import BeautifulSoup
import re

# 🔹 1. Définition des primers et de l'espèce
primer_fwd = "GGCAACACTCTCGGAGACAA"
primer_rev = "GGAAAGATCCCAGCAGCAGT"
species = "Homo+sapiens"

# 🔹 2. Construction de l'URL de requête
url = (
    "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?"
    "CMD=request&CON_ANEAL_OLIGO=50.0&CON_DNTPS=0.6&DIVA_CATIONS=1.5&EVALUE=30000"
    "&GC_CLAMP=0&HITSIZE=50000&LOW_COMPLEXITY_FILTER=on&MAX_CANDIDATE_PRIMER=500"
    "&MAX_INTRON_SIZE=1000000&MAX_TARGET_PER_TEMPLATE=100&MAX_TARGET_SIZE=4000"
    "&MIN_INTRON_SIZE=1000&MISMATCH_REGION_LENGTH=5&MONO_CATIONS=50.0"
    f"&ORGANISM={species}&OVERLAP_3END=4&OVERLAP_5END=7&POLYX=5"
    "&PRIMER_3END_SPECIFICITY_MISMATCH=1&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0"
    "&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0"
    "&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18"
    "&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50"
    "&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0"
    f"&PRIMER_LEFT_INPUT={primer_fwd}&PRIMER_MAX_DIFF_TM=3&PRIMER_MAX_END_GC=5"
    "&PRIMER_MAX_END_STABILITY=9&PRIMER_MAX_GC=80.0&PRIMER_MAX_HAIRPIN_TH=24.0"
    "&PRIMER_MAX_SELF_ANY_TH=45.0&PRIMER_MAX_SELF_END_TH=35.0"
    "&PRIMER_MAX_SIZE=25&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00"
    "&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00&PRIMER_MAX_TM=63.0"
    "&PRIMER_MIN_GC=20.0&PRIMER_MIN_SIZE=15&PRIMER_MIN_TM=57.0"
    "&PRIMER_MISPRIMING_LIBRARY=AUTO&PRIMER_NUM_RETURN=10&PRIMER_ON_SPLICE_SITE=0"
    "&PRIMER_OPT_SIZE=20&PRIMER_OPT_TM=60.0&PRIMER_PAIR_MAX_COMPL_ANY=8.00"
    "&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0&PRIMER_PAIR_MAX_COMPL_END=3.00"
    "&PRIMER_PAIR_MAX_COMPL_END_TH=35.0&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00"
    "&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00&PRIMER_PRODUCT_MAX=1000"
    "&PRIMER_PRODUCT_MIN=70"
    f"&PRIMER_RIGHT_INPUT={primer_rev}&PRIMER_SPECIFICITY_DATABASE=refseq_mrna"
    "&SALT_FORMULAR=1&SEARCH_SPECIFIC_PRIMER=on&SEARCHMODE=0&SELF_ANY=8.00"
    "&SELF_END=3.00&SHOW_SVIEWER=on&SPLICE_SITE_OVERLAP_3END=4"
    "&SPLICE_SITE_OVERLAP_3END_MAX=8&SPLICE_SITE_OVERLAP_5END=7&TM_METHOD=1"
    "&TOTAL_MISMATCH_IGNORE=6&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1&UNGAPPED_BLAST=on"
    "&USER_TYPE=2&WORD_SIZE=7"
)

# 🔹 3. Envoi de la requête POST et récupération de la réponse HTML
session = requests.Session()
response = session.get(url)

# Vérification si la requête a réussi
if response.status_code != 200:
    print(f"❌ Erreur lors de la requête : {response.status_code}")
    exit()

# 🔹 4. Extraction du JOB ID depuis le HTML
soup = BeautifulSoup(response.text, "html.parser")
# Méthode 1 : Rechercher un champ input caché avec le JOB ID
job_id_input = soup.find("input", {"name": "job_key"})
if job_id_input:
    job_key = job_id_input.get("value")
else:
    # Méthode 2 : Chercher le JOB ID dans du JavaScript
    match = re.search(r"job_key=([A-Za-z0-9_-]+)", response.text)
    if match:
        job_key = match.group(1)
    else:
        print("❌ Impossible de récupérer le JOB ID.")
        exit()

print(f"✅ JOB ID récupéré : {job_key}")

# 🔹 5. Vérification périodique de l'état du job
status_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key={job_key}&CMD=get"

while True:
    time.sleep(10)  # Attendre 10 secondes avant chaque vérification
    status_response = session.get(status_url)

    soup = BeautifulSoup(status_response.text, "html.parser")
    print(soup)

    # Vérifier si les résultats sont disponibles
    primer_pair = soup.find("a", {"name": "0"})
    forward_primer = soup.find("th", string="Forward primer")
    reverse_primer = soup.find("th", string="Reverse primer")

    if primer_pair and forward_primer and reverse_primer:
        print("✅ Résultats prêts !")
        break
    else:
        print("⏳ En attente des résultats...")

# 🔹 6. Extraction des résultats PCR in-silico
results = []
amplicon_sizes = soup.find_all("td", string=re.compile(r"\d+\s*bp"))

for row in amplicon_sizes:
    results.append(row.text.strip())

if results:
    print(f"✅ Amplicons trouvés : {results}")
else:
    print("❌ Aucun amplicon trouvé.")
