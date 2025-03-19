import requests
import re

BASE_URL = "https://enzymefinder.neb.com"

try:
    # Récupérer la page principale
    response = requests.get(BASE_URL, timeout=10)
    response.raise_for_status()

    # Chercher le script correspondant à "main-XXXXXX.js"
    match = re.search(r'src="(/scripts/main-[a-f0-9]+\.js)"', response.text)

    if match:
        script_url = BASE_URL + match.group(1)
        print(f"Tentative de téléchargement : {script_url}")

        # Télécharger le fichier JS
        script_response = requests.get(script_url, timeout=10)
        script_response.raise_for_status()

        script_content = script_response.text
        print("Script téléchargé avec succès.")
    else:
        print("⚠️ Aucun script correspondant trouvé !")
        script_content = None

except requests.RequestException as e:
    print(f"Erreur lors de la récupération : {e}")
    script_content = None
