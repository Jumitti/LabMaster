import re

import requests
import streamlit as st
import pandas as pd

from pages.enzymes_utils import NEB


@st.cache_resource
def get_PROMEGA_enzyme_data():
    url = 'https://scv10mr-cdnpre-p-cus-00.azureedge.net/-/media/files/javascript/retool/enzymes.js'
    local_file_path = 'pages/enzymes_utils/PROMEGA_enz_db.js'

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        script_content = response.text
    except (requests.RequestException, requests.Timeout):
        with open(local_file_path, 'r') as file:
            script_content = file.read()
        st.warning("⚠️ The Promega database seems inaccessible. A local version has been loaded and may have differences with the online version")

    enzyme_pattern = re.compile(r'\{(.*?)},', re.DOTALL)
    enzyme_blocks = enzyme_pattern.findall(script_content)
    enzymes = []
    for block in enzyme_blocks:
        name_match = re.search(r"name:'(.*?)'", block)
        temperature_match = re.search(r"temperature:'(.*?)'", block)
        star_match = re.search(r'star:\[(.*?)]', block)
        bufferActivity_match = re.search(r"bufferActivity:\{(.*?)}", block, re.DOTALL)

        name = name_match.group(1) if name_match else ""
        temperature = temperature_match.group(1) if temperature_match else ""
        star = [s.strip().strip("'") for s in star_match.group(1).split(',') if
                s.strip().strip("'")] if star_match else []

        bufferActivity = {}
        if bufferActivity_match:
            bufferActivity_str = bufferActivity_match.group(1)
            bufferActivity_items = bufferActivity_str.split(',')
            for item in bufferActivity_items:
                if ':' in item:
                    buffer, activity = item.split(':')
                    bufferActivity[buffer.strip()] = int(activity.strip())

        NEB_enzymes = NEB.get_NEB_enzyme_data()
        NEB_dict = {enzyme['Enzyme']: enzyme for enzyme in NEB_enzymes}
        heatInactivation_PROMEGA_csv = pd.read_csv('pages/enzymes_utils/PROMEGA_enz_hi_db.csv', sep=';')
        heatInactivation_PROMEGA_dict = dict(zip(heatInactivation_PROMEGA_csv['Enzyme'], heatInactivation_PROMEGA_csv['HeatInactivation']))

        heat_inactivation_symbols = {
            "yes": "✅",
            "partial": "⚠️",
            "no": "❌"
        }
        heatInactivation = heat_inactivation_symbols.get(heatInactivation_PROMEGA_dict.get(name, ""), "")

        enzyme = {
            'Enzyme': name,
            'Type': f"{NEB_dict[name]['Type']}" if name in NEB_dict else "",
            'Subtype': f"{NEB_dict[name]['Subtype']}" if name in NEB_dict else "",
            'Enzyme Type': f"{NEB_dict[name]['Enzyme Type']}" if name in NEB_dict else "",
            'Overhang': f"{NEB_dict[name]['Overhang']}" if name in NEB_dict else "",
            'Sequence': f"{NEB_dict[name]['Sequence']}" if name in NEB_dict else "",
            'Overhang sequence': f"{NEB_dict[name]['Overhang sequence']}" if name in NEB_dict else "",
            'Star': star,
            'Buffer A': f"{bufferActivity.get('A', '')}*" if "A" in star else bufferActivity.get('A', ''),
            'Buffer B': f"{bufferActivity.get('B', '')}*" if "B" in star else bufferActivity.get('B', ''),
            'Buffer C': f"{bufferActivity.get('C', '')}*" if "C" in star else bufferActivity.get('C', ''),
            'Buffer D': f"{bufferActivity.get('D', '')}*" if "D" in star else bufferActivity.get('D', ''),
            'Buffer E': f"{bufferActivity.get('E', '')}*" if "E" in star else bufferActivity.get('E', ''),
            'Buffer F': f"{bufferActivity.get('F', '')}*" if "E" in star else bufferActivity.get('F', ''),
            'Buffer G': f"{bufferActivity.get('G', '')}*" if "H" in star else bufferActivity.get('G', ''),
            'Buffer H': f"{bufferActivity.get('H', '')}*" if "H" in star else bufferActivity.get('H', ''),
            'Buffer J': f"{bufferActivity.get('J', '')}*" if "H" in star else bufferActivity.get('J', ''),
            'Buffer K': f"{bufferActivity.get('K', '')}*" if "H" in star else bufferActivity.get('K', ''),
            'Buffer MultiCore': f"{bufferActivity.get('MultiCore', '')}*" if "MultiCore" in star else bufferActivity.get('MultiCore', ''),
            'Incubation Temp (°C)': temperature,
            'Heat Inactivation': heatInactivation,
            'Multiple sites': f"{NEB_dict[name]['Multiple sites']}" if name in NEB_dict else "",
            'CpG': f"{NEB_dict[name]['CpG']}" if name in NEB_dict else "",
            'Dam': f"{NEB_dict[name]['Dam']}" if name in NEB_dict else "",
            'Dcm': f"{NEB_dict[name]['Dcm']}" if name in NEB_dict else "",
            'Methylation dependence': NEB_dict[name]['Methylation dependence'] if name in NEB_dict else "",
            'Methylation sensitivity': NEB_dict[name]['Methylation sensitivity'] if name in NEB_dict else [],
            'Isoschizomers': NEB_dict[name]['Isoschizomers'] if name in NEB_dict else [],
            'url': f'https://www.promega.com/products/cloning-and-dna-markers/restriction-enzymes/{name.lower()}'
        }

        enzymes.append(enzyme)

    return enzymes