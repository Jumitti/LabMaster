import re

import requests
import streamlit as st


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
def get_NEB_enzyme_data():
    url = 'https://enzymefinder.neb.com/scripts/main-b4fe94d919.js'
    local_file_path = 'pages/enzymes_utils/NEB_enz_db.js'

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        script_content = response.text
    except (requests.RequestException, requests.Timeout):
        with open(local_file_path, 'r') as file:
            script_content = file.read()
        st.warning(
            "⚠️ The NEB database seems inaccessible. A local version has been loaded and may have differences with the online version")

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
        star1_match = re.search(r'star1:(true|false|!0|!1)', block)
        star2_match = re.search(r'star2:(true|false|!0|!1)', block)
        star3_match = re.search(r'star3:(true|false|!0|!1)', block)
        star4_match = re.search(r'star4:(true|false|!0|!1)', block)
        cpg_match = re.search(r'cpg:(true|false|!0|!1)', block)
        dam_match = re.search(r'dam:(true|false|!0|!1)', block)
        dcm_match = re.search(r'dcm:(true|false|!0|!1)', block)
        heatInactivationTemp_match = re.search(r'heatInactivationTemp:(\d+)', block)
        heatInactivationTime_match = re.search(r'heatInactivationTime:(\d+)', block)
        incubateTemp_match = re.search(r'incubateTemp:(\d+)', block)
        timeSaver_match = re.search(r'timeSaver:(true|false|!0|!1)', block)
        url_match = re.search(r'url:"(.*?)"', block)
        mul_match = re.search(r'mul:(true|false|!0|!1)', block)
        atp_match = re.search(r'supplement:{.*?atp:(\d+)', block)
        bsa_match = re.search(r'supplement:{.*?bsa:(\d+)', block)
        dtt_match = re.search(r'supplement:{.*?dtt:(\d+)', block)
        enzact_match = re.search(r'supplement:{.*?enzact:(\d+)', block)
        hfEnzyme_match = re.search(r'hfEnzyme:(true|false|!0|!1)', block)

        displayName = displayName_match.group(1) if displayName_match else ""
        overhangType = overhangType_match.group(1) if overhangType_match else ""
        HFVersionAvailable = HFVersionAvailable_match.group(1) if HFVersionAvailable_match else ""
        soldByNEB = soldByNEB_match.group(1) if soldByNEB_match else ""
        subtype = subtype_match.group(1) if subtype_match else ""
        type_ = type_match.group(1) if type_match else ""
        md = md_match.group(1) if md_match else ""
        ms = [ms_enz.strip().strip('"') for ms_enz in ms_match.group(1).split(',') if ms_enz.strip().strip('"')] if ms_match else [""]
        compactFormattedSeq = convert_to_subscript_superscript(
            compactFormattedSeq_match.group(1)) if compactFormattedSeq_match else ""
        isoschizomers = [iso.strip().strip('"') for iso in
                         isoschizomers_match.group(1).split(',') if iso.strip().strip('"')] if isoschizomers_match else [""]
        nicking = nicking_match.group(1) if nicking_match else ""
        overhangSeq = overhangSeq_match.group(1) if overhangSeq_match else ""
        buf1 = buf1_match.group(1) if buf1_match else ""
        buf2 = buf2_match.group(1) if buf2_match else ""
        buf3 = buf3_match.group(1) if buf3_match else ""
        buf4 = buf4_match.group(1) if buf4_match else ""
        star1 = star1_match.group(1) if star1_match else ""
        star2 = star2_match.group(1) if star2_match else ""
        star3 = star3_match.group(1) if star3_match else ""
        star4 = star4_match.group(1) if star4_match else ""
        cpg = cpg_match.group(1) if cpg_match else ""
        dam = dam_match.group(1) if dam_match else ""
        dcm = dcm_match.group(1) if dcm_match else ""
        heatInactivationTemp = heatInactivationTemp_match.group(1) if heatInactivationTemp_match else ""
        heatInactivationTime = heatInactivationTime_match.group(1) if heatInactivationTime_match else ""
        incubateTemp = incubateTemp_match.group(1) if incubateTemp_match else ""
        timeSaver = timeSaver_match.group(1) if timeSaver_match else ""
        url = url_match.group(1) if url_match else ""
        mul = mul_match.group(1) if mul_match else ""
        atp = atp_match.group(1) if atp_match else ""
        bsa = bsa_match.group(1) if bsa_match else ""
        dtt = dtt_match.group(1) if dtt_match else ""
        enzact = enzact_match.group(1) if enzact_match else ""
        hfEnzyme = hfEnzyme_match.group(1) if hfEnzyme_match else ""

        if HFVersionAvailable != "":
            HFVersionAvailable = "✅" if HFVersionAvailable in ["!0", "true"] else "❌"
        if soldByNEB != "":
            soldByNEB = "✅" if soldByNEB in ["!0", "true"] else "❌"
        if nicking != "":
            if nicking in ["!0", "true"]:
                nicking = "✅"
                subtype = "Nicking"
            else:
                nicking = "❌"
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
        if hfEnzyme != "":
            hfEnzyme = "Yes" if hfEnzyme in ["!0", "true"] else ""

        stars = [star1, star2, star3, star4]
        for i in range(len(stars)):
            if stars[i] != "":
                stars[i] = "*" if stars[i] in ["!0", "true"] else ""
        star1, star2, star3, star4 = stars

        ms = [item.replace("-", "untested") for item in ms]

        enzyme = {
            'Enzyme': displayName,
            'Type': type_,
            'Subtype': subtype,
            'Enzyme Type': f'{type_} {subtype}',
            'Overhang': overhangType,
            'Sequence': compactFormattedSeq,
            'Overhang sequence': overhangSeq,
            'NEBuffer r1.1': f"{buf1}{star1}" if buf1 != "" else "",
            'NEBuffer r2.1': f"{buf2}{star2}" if buf2 != "" else "",
            'NEBuffer r3.1': f"{buf3}{star3}" if buf3 != "" else "",
            'rCutSmart': f"{buf4}{star4}" if buf4 != "" else "",
            'ATP': atp,
            'BSA': bsa,
            'DTT': dtt,
            'ENZACT': enzact,
            'Supplement': f'{"ATP" if atp not in ["0", ""] else "" "BSA" if bsa not in ["0", ""] else "" "DTT" if dtt not in ["0", ""] else "" "Enzyme Activator" if enzact not in ["0", ""] else ""}',
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
            'HF version available': HFVersionAvailable,
            'HF version': hfEnzyme,
            'Sold by NEB': soldByNEB,
            'url': url,

        }

        enzymes.append(enzyme)

    return enzymes