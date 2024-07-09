import re
import pandas as pd
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
        buf1_match = re.search(r'buf1:(\d+)', block)
        buf2_match = re.search(r'buf2:(\d+)', block)
        buf3_match = re.search(r'buf3:(\d+)', block)
        buf4_match = re.search(r'buf4:(\d+)', block)
        cpg_match = re.search(r'cpg:(true|false|!0|!1)', block)
        dam_match = re.search(r'dam:(true|false|!0|!1)', block)
        dcm_match = re.search(r'dcm:(true|false|!0|!1)', block)
        heatInactivationTemp_match = re.search(r'heatInactivationTemp:(\d+)', block)
        heatInactivationTime_match = re.search(r'heatInactivationTime:(\d+)', block)
        incubateTemp_match = re.search(r'incubateTemp:(\d+)', block)
        timeSaver_match = re.search(r'timeSaver:(true|false|!0|!1)', block)
        url_match = re.search(r'url:"(.*?)"', block)
        mul_match = re.search(r'mul:(true|false|!0|!1)', block)

        displayName = displayName_match.group(1) if displayName_match else ""
        overhangType = overhangType_match.group(1) if overhangType_match else ""
        HFVersionAvailable = HFVersionAvailable_match.group(1) if HFVersionAvailable_match else ""
        soldByNEB = soldByNEB_match.group(1) if soldByNEB_match else ""
        subtype = subtype_match.group(1) if subtype_match else ""
        type_ = type_match.group(1) if type_match else ""
        md = md_match.group(1) if md_match else ""
        ms = [ms_enz.strip().strip('"') for ms_enz in ms_match.group(1).split(',')] if ms_match else [""]
        compactFormattedSeq = convert_to_subscript_superscript(compactFormattedSeq_match.group(1)) if compactFormattedSeq_match else ""
        isoschizomers = [iso.strip().strip('"') for iso in isoschizomers_match.group(1).split(',')] if isoschizomers_match else [""]
        nicking = nicking_match.group(1) if nicking_match else ""
        overhangSeq = overhangSeq_match.group(1) if overhangSeq_match else ""
        buf1 = buf1_match.group(1) if buf1_match else ""
        buf2 = buf2_match.group(1) if buf2_match else ""
        buf3 = buf3_match.group(1) if buf3_match else ""
        buf4 = buf4_match.group(1) if buf4_match else ""
        cpg = cpg_match.group(1) if cpg_match else ""
        dam = dam_match.group(1) if dam_match else ""
        dcm = dcm_match.group(1) if dcm_match else ""
        heatInactivationTemp = heatInactivationTemp_match.group(1) if heatInactivationTemp_match else ""
        heatInactivationTime = heatInactivationTime_match.group(1) if heatInactivationTime_match else ""
        incubateTemp = incubateTemp_match.group(1) if incubateTemp_match else ""
        timeSaver = timeSaver_match.group(1) if timeSaver_match else ""
        url = url_match.group(1) if url_match else ""
        mul = mul_match.group(1) if mul_match else ""

        if HFVersionAvailable != "":
            HFVersionAvailable = "✅" if HFVersionAvailable in ["!0", "true"] else "❌"
        if soldByNEB != "":
            soldByNEB = "✅" if soldByNEB in ["!0", "true"] else "❌"
        if nicking != "":
            nicking = "✅" if nicking in ["!0", "true"] else "❌"
        if overhangSeq == "":
            overhangSeq = overhangType
        elif overhangType == "5-Prime":
            overhangSeq = f"5' {overhangSeq}"
        elif overhangType == "3-Prime":
            overhangSeq = f"3' {overhangSeq}"
        if cpg != "":
            cpg = "✅" if cpg in ["!0", "true"] else "❌"
        if dam != "":
            dam = "✅" if dam in ["!0", "true"] else "❌"
        if dcm != "":
            dcm = "✅" if dcm in ["!0", "true"] else "❌"
        if timeSaver != "":
            timeSaver = "✅" if timeSaver in ["!0", "true"] else "❌"
        if mul != "":
            mul = "requires 2 sites" if mul in ["!0", "true"] else ""

        ms = [item.replace("-", "untested") for item in ms]

        enzyme = {
            'Enzyme': displayName,
            'Type': type_,
            'Subtype': subtype,
            'Enzyme Type': f'{type_} {subtype}',
            'Overhang': overhangType,
            'Sequence': compactFormattedSeq,
            'Overhang sequence': overhangSeq,
            'NEBuffer r1.1': buf1,
            'NEBuffer r2.1': buf2,
            'NEBuffer r3.1': buf3,
            'rCutSmart': buf4,
            'Incubation Temp (°C)': incubateTemp,
            'TimeSaver': timeSaver,
            'heatInactivationTemp': heatInactivationTemp,
            'heatInactivationTime': heatInactivationTime,
            'Heat Inactivation': f'{heatInactivationTemp}°C for {heatInactivationTime}min' if heatInactivationTemp not in [
                "", "0"] and heatInactivationTime not in ["", "0"] else "",
            'Multiple sites': mul,
            'CpG': cpg,
            'Dam': dam,
            'Dcm': dcm,
            'Methylation dependence': md,
            'Methylation sensitivity': ms,
            'Isoschizomers': isoschizomers,
            'Nicking': nicking,
            'HF version': HFVersionAvailable,
            'Sold by NEB': soldByNEB,
            'url': url

        }

        enzymes.append(enzyme)

    return enzymes


# Page config
page_config()

st.warning('⚠️ This is a DEMO still in testing... Enzymes are from NEB database')

enzymes = get_enzyme_data()

enzymes_df = pd.DataFrame(enzymes)

columns_to_exclude = ['Type', 'Subtype', 'heatInactivationTemp', 'heatInactivationTime']
enzymes_df_filtered = enzymes_df.drop(columns=columns_to_exclude)

NEB_only = st.toggle("Show NEB Enzymes Only", value=True)
if NEB_only:
    enzymes_df_filtered = enzymes_df_filtered[enzymes_df['Sold by NEB'] == "✅"]

type_only = st.selectbox("Enzyme type", ["Any", "Type II", "Homing"])
if type_only != "Any":
    enzymes_df_filtered = enzymes_df_filtered[enzymes_df['Type'].str.contains(type_only, case=False, na=False)]

st.write(f"{len(enzymes_df_filtered)} enzymes available")
st.dataframe(enzymes_df_filtered,
             column_config={
                 "url": st.column_config.LinkColumn("URL")}, hide_index=True)

selected_names = st.multiselect("Enzymes", [enzyme['Enzyme'] for enzyme in enzymes])

selected_enzymes = [enzyme for enzyme in enzymes if enzyme['Enzyme'] in selected_names]

# Display the selected enzymes
for enzyme in selected_enzymes:
    st.write(enzyme)

