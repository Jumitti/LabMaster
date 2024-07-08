import re

import requests
import streamlit as st

from utils.page_config import page_config


@st.cache
def get_enzyme_data():
    url = 'https://enzymefinder.neb.com/scripts/main-b4fe94d919.js'

    response = requests.get(url)
    response.raise_for_status()

    script_content = response.text

    pattern = re.compile(r'displayName:"(.*?)".*?overhangType:"(.*?)".*?soldByNEB:(true|false|!1|!0)', re.DOTALL)
    matches = pattern.findall(script_content)

    enzymes = []
    for match in matches:
        displayName, overhangType, soldByNEB = match
        soldByNEB = "Yes" if soldByNEB in ["!1", "true"] else "No"
        enzymes.append({'displayName': displayName, 'overhangType': overhangType, 'soldByNEB': soldByNEB})

    return enzymes


# Page config
page_config()

enz = get_enzyme_data()

if enz:
    st.write(enz)
else:
    st.write('No data')
