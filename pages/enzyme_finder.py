import re

import requests
import streamlit as st

from utils.page_config import page_config


def convert_to_subscript_superscript(seq):
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    superscript_map = str.maketrans("mh", "ᵐʰ")
    seq = seq.replace("^", "↓")
    seq = seq.replace("_", "↑")
    seq = seq.replace("n", "N")
    seq = seq.translate(subscript_map)
    seq = seq.translate(superscript_map)
    return seq

@st.cache_resource
def get_enzyme_data():
    url = 'https://enzymefinder.neb.com/scripts/main-b4fe94d919.js'

    response = requests.get(url)
    response.raise_for_status()

    script_content = response.text

    # Regex pattern to capture each enzyme block
    enzyme_pattern = re.compile(r'\{[^{}]*displayName:".*?"[^{}]*}', re.DOTALL)
    enzyme_blocks = enzyme_pattern.findall(script_content)

    enzymes = []
    for block in enzyme_blocks:
        displayName_match = re.search(r'displayName:"(.*?)"', block)
        overhangType_match = re.search(r'overhangType:"(.*?)"', block)
        HFVersionAvailable_match = re.search(r'HFVersionAvailable:(true|false|!0|!1)', block)
        soldByNEB_match = re.search(r'soldByNEB:(true|false|!0|!1)', block)
        subtype_match = re.search(r'subtype:"(.*?)"', block)
        type_match = re.search(r'(?<!sub)type:"(.*?)"', block)
        md_match = re.search(r'md:"(.*?)"', block)
        ms_match = re.search(r'ms:\[(.*?)]', block)
        compactFormattedSeq_match = re.search(r'compactFormattedSeq:"(.*?)"', block)
        isoschizomers_match = re.search(r'isoschizomers:\[(.*?)]', block)
        nicking_match = re.search(r'nicking:(true|false|!0|!1)', block)
        overhangSeq_match = re.search(r'overhangSeq:"(.*?)"', block)

        displayName = displayName_match.group(1) if displayName_match else "N/A"
        overhangType = overhangType_match.group(1) if overhangType_match else "N/A"
        HFVersionAvailable = HFVersionAvailable_match.group(1) if HFVersionAvailable_match else "N/A"
        soldByNEB = soldByNEB_match.group(1) if soldByNEB_match else "N/A"
        subtype = subtype_match.group(1) if subtype_match else "N/A"
        type_ = type_match.group(1) if type_match else "N/A"
        md = md_match.group(1) if md_match else "N/A"
        ms = ms_match.group(1).split(',') if ms_match else []
        compactFormattedSeq = convert_to_subscript_superscript(compactFormattedSeq_match.group(1)) if compactFormattedSeq_match else "N/A"
        isoschizomers = isoschizomers_match.group(1).split(',') if isoschizomers_match else []
        nicking = nicking_match.group(1) if nicking_match else "N/A"
        overhangSeq = overhangSeq_match.group(1) if overhangSeq_match else "N/A"

        if HFVersionAvailable != "N/A":
            HFVersionAvailable = "Yes" if HFVersionAvailable in ["!0", "true"] else "No"
        if soldByNEB != "N/A":
            soldByNEB = "Yes" if soldByNEB in ["!0", "true"] else "No"
        if nicking != "N/A":
            nicking = "Yes" if nicking in ["!0", "true"] else "No"
        if overhangSeq != "N/A":
            if overhangSeq == "":
                overhangSeq = overhangType
            elif overhangType == "5-Prime":
                overhangSeq = f"5' {overhangSeq}"
            elif overhangType == "3-Prime":
                overhangSeq = f"3' {overhangSeq}"

        ms = [item.replace("-", "untested") for item in ms]

        enzyme = {
            'displayName': displayName,
            'overhangType': overhangType,
            'HFVersionAvailable': HFVersionAvailable,
            'soldByNEB': soldByNEB,
            'subtype': subtype,
            'type': type_,
            'md': md,
            'ms': ms,
            'compactFormattedSeq': compactFormattedSeq,
            'isoschizomers':  isoschizomers,
            'nicking': nicking,
            'overhangSeq': overhangSeq

        }

        enzymes.append(enzyme)

    return enzymes


# Page config
page_config()

enzymes = get_enzyme_data()

st.write(f"{len(enzymes)} enzymes available")
selected_names = st.multiselect("Enzymes", [enzyme['displayName'] for enzyme in enzymes])

selected_enzymes = [enzyme for enzyme in enzymes if enzyme['displayName'] in selected_names]

st.write(enzymes)
# Display the selected enzymes
for enzyme in selected_enzymes:
    st.write(enzyme)

