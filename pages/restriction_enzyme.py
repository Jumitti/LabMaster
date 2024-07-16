import re
import pandas as pd
import requests
import streamlit as st

from utils.page_config import page_config
from pages.enzymes_utils import NEB, PROMEGA


def display_table(db="NEB"):
    if db == "NEB":
        enzymes = NEB.get_NEB_enzyme_data()
    elif db == 'PROMEGA':
        enzymes = PROMEGA.get_PROMEGA_enzyme_data()
    enzymes_df = pd.DataFrame(enzymes)

    col1, col2, col3 = st.columns(3, gap="small")
    overhang_only = col1.selectbox("Overhang type", ["Any", "Blunt", "5-Prime", '3-Prime'], key=f"{db}_overhang")
    if overhang_only != "Any":
        enzymes_df = enzymes_df[enzymes_df['Overhang'].str.contains(overhang_only, case=False, na=False)]

    type_only = col2.selectbox("Enzyme type", ["Any", "Type II", "Homing", 'Nicking'], key=f"{db}_enzyme_type")
    if type_only != "Any":
        if type_only != "Nicking":
            enzymes_df = enzymes_df[enzymes_df['Type'].str.contains(type_only, case=False, na=False)]
        elif type_only == "Nicking":
            enzymes_df = enzymes_df[enzymes_df['Subtype'].str.contains(type_only, case=False, na=False)]

    subtype_only = col3.selectbox("Enzyme subtype",
                                  ["Any", "A", "a", "B", "c", "E", "F", "G", "H", "M", "P", "S", "unknown"],
                                  key=f"{db}_enzyme_subtype")
    if subtype_only != "Any":
        enzymes_df = enzymes_df[enzymes_df['Subtype'].str.contains(subtype_only, case=False, na=False)]

    if db == "NEB":
        enzymes_df = enzymes_df[enzymes_df['Sold by NEB'] == "✅"]
        HF_only = st.toggle("Show HF NEB Enzymes Only", value=False, key="HF",
                            help="For over 45 years, New England Biolabs has been developing innovative solutions for molecular biology applications. The respected leader in the field of restriction enzyme biology, NEB has developed a line of High-Fidelity (HF®) Restriction Enzymes. These engineered enzymes have the same specificity as the native enzyme, with the added benefit of reduced star activity, rapid digestion (5-15 minutes) and 100% activity in rCutSmart™ and CutSmart® Buffer.\n\nEngineered with performance in mind, HF restriction enzymes are fully active under a broader range of conditions, minimizing off-target products, while offering flexibility in experimental design.")
        if HF_only:
            enzymes_df = enzymes_df[enzymes_df['HF version'] == "Yes"]
        columns_to_exclude = ['Type', "Subtype", 'heatInactivationTemp', 'heatInactivationTime', 'ATP', 'BSA', 'DTT',
                              'ENZACT', 'HF version', 'Nicking', 'Sold by NEB']
    elif db == "PROMEGA":
        enzymes_df = enzymes_df[enzymes_df['Buffer A'] != ""]
        columns_to_exclude = ['Type', "Subtype", "Star"]

    enzymes_df_filtered = enzymes_df.drop(columns=columns_to_exclude)

    st.dataframe(enzymes_df_filtered.sort_values(by='Enzyme').style.map(
        highlight_buffer, subset=['NEBuffer r1.1', 'NEBuffer r2.1', 'NEBuffer r3.1', 'rCutSmart'] if
        db == "NEB" else ['Buffer A', 'Buffer B', 'Buffer C', 'Buffer D', 'Buffer E', 'Buffer F', 'Buffer G',
                          'Buffer H', 'Buffer J', 'Buffer K', 'Buffer MultiCore']), column_config={"url": st.column_config.LinkColumn("URL")}, hide_index=True, key=f"{db}_dataframe")
    st.write(
        f"**{len(enzymes_df_filtered)}** enzymes available | Buffers are in percentage of **effectiveness** - if there is **a *(star)**, the enzyme has **star** activity in this buffer. {' | NEB has **TimeSaver** enzymes which digest in **5-15 minutes**.' if db == 'NEB' else ''}")

    return enzymes


def select_enzyme(enzymes, db="NEB"):
    st.divider()
    selected_names = st.multiselect("Select enzymes", [enzyme['Enzyme'] for enzyme in enzymes])
    if selected_names:
        selected_enzymes = [enzyme for enzyme in enzymes if enzyme['Enzyme'] in selected_names]
        enzymes_df = pd.DataFrame(selected_enzymes)

        if db == "NEB":
            columns_to_exclude = ['Type', "Subtype", 'heatInactivationTemp', 'heatInactivationTime', 'ATP', 'BSA',
                                  'DTT',
                                  'ENZACT', 'HF version', 'Nicking', 'Sold by NEB']
        elif db == "PROMEGA":
            columns_to_exclude = ['Type', "Subtype", "Star"]
        enzymes_df_filtered = enzymes_df.drop(columns=columns_to_exclude)

        st.dataframe(enzymes_df_filtered.sort_values(by='Enzyme').style.map(
            highlight_buffer, subset=['NEBuffer r1.1', 'NEBuffer r2.1', 'NEBuffer r3.1', 'rCutSmart'] if
            db == "NEB" else ['Buffer A', 'Buffer B', 'Buffer C', 'Buffer D', 'Buffer E', 'Buffer F', 'Buffer G',
                              'Buffer H',
                              'Buffer J', 'Buffer K', 'Buffer MultiCore']), column_config={
            "url": st.column_config.LinkColumn("URL")}, hide_index=True, key=f"{db}_dataframe")

        return selected_enzymes


def digestion_protocols(enzymes, db="NEB"):
    st.divider()
    col1dp, col2dp, col3dp = st.columns(3, gap="small")
    dna_reaction = col1dp.number_input("DNA per reaction (ng)", value=1000.00, step=0.01, min_value=0.00,
                                       key=f"{db}_dna_reaction")
    dna_conc = col2dp.number_input("DNA concentration (ng/µL)", value=100.00, step=0.01, min_value=0.00,
                                   key=f"{db}_dna_concentration")
    samples = col3dp.number_input("DNA per reaction (ng)", value=2, step=1, min_value=1, key=f"{db}_samples")

    if db == "NEB":
        protocols_output = [
                               {"Component": f"DNA ({dna_reaction / 1000}µg)",
                                "For 1 reaction (µL)": dna_reaction / dna_conc,
                                f"For {samples} reactions (µL)": samples * dna_reaction / dna_conc},
                               {"Component": "Buffer (µL)", "For 1 reaction (µL)": 5,
                                f"For {samples} reactions (µL)": samples * 5}
                           ] + [
                               {"Component": f"{enzyme['Enzyme']} (µL)", "For 1 reaction (µL)": 1,
                                f"For {samples} reactions (µL)": samples * 1,
                                'NEBuffer r1.1': enzyme['NEBuffer r1.1'], 'NEBuffer r2.1': enzyme['NEBuffer r2.1'],
                                'NEBuffer r3.1': enzyme['NEBuffer r3.1'],
                                'rCutSmart': enzyme['rCutSmart'], 'TimeSaver': enzyme['TimeSaver']} for enzyme in
                               enzymes
                           ] + [
                               {"Component": f"H2O (μL)",
                                "For 1 reaction (µL)": 50 - len(enzymes) - 5 - dna_reaction / dna_conc,
                                f"For {samples} reactions (µL)": samples * (
                                        50 - samples * len(enzymes) - 5 - dna_reaction / dna_conc)},
                               {"Component": "Total (µL)", "For 1 reaction (µL)": 50 if 50 - len(
                                   enzymes) - 5 - dna_reaction / dna_conc >= 0 else 50 - (
                                       50 - len(enzymes) - 5 - dna_reaction / dna_conc),
                                f"For {samples} reactions (µL)": samples * 50 if 50 - len(
                                    enzymes) - 5 - dna_reaction / dna_conc >= 0 else samples * (
                                        50 - (50 - len(enzymes) - 5 - dna_reaction / dna_conc))}]
    if db == "PROMEGA":
        vol_enz = 1 if len(enzymes) > 1 else 2
        protocols_output = [
                               {"Component": f"DNA ({dna_reaction / 1000}µg)",
                                "For 1 reaction (µL)": dna_reaction / dna_conc,
                                f"For {samples} reactions (µL)": samples * dna_reaction / dna_conc},
                               {"Component": "Buffer (µL)", "For 1 reaction (µL)": 2,
                                f"For {samples} reactions (µL)": samples * 2},
                               {"Component": "BSA 10X (µL)", "For 1 reaction (µL)": 2,
                                f"For {samples} reactions (µL)": samples * 2}
                           ] + [
                               {"Component": f"{enzyme['Enzyme']} (µL)", "For 1 reaction (µL)": vol_enz,
                                f"For {samples} reactions (µL)": samples * vol_enz,
                                'Buffer A': f"{enzyme['Buffer A']}", 'Buffer B': f"{enzyme['Buffer B']}",
                                'Buffer C': f"{enzyme['Buffer C']}",
                                'Buffer D': f"{enzyme['Buffer D']}", 'Buffer E': f"{enzyme['Buffer E']}",
                                'Buffer F': f"{enzyme['Buffer F']}",
                                'Buffer G': f"{enzyme['Buffer G']}", 'Buffer H': f"{enzyme['Buffer H']}",
                                'Buffer J': f"{enzyme['Buffer J']}",
                                'Buffer K': f"{enzyme['Buffer K']}",
                                'Buffer MultiCore': f"{enzyme['Buffer MultiCore']}",
                                } for enzyme in enzymes
                           ] + [
                               {"Component": f"H2O (μL)",
                                "For 1 reaction (µL)": 20 - len(enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc,
                                f"For {samples} reactions (µL)": samples * (
                                        20 - samples * len(enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc)},
                               {"Component": "Total (µL)",
                                "For 1 reaction (µL)": 20 if 20 - len(
                                    enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc >= 0 else 20 - (
                                        20 - len(enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc),
                                f"For {samples} reactions (µL)": samples * 20 if 20 - len(
                                    enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc >= 0 else samples * (20 - (
                                        20 - len(enzymes) * vol_enz - 2 - 2 - dna_reaction / dna_conc))}]

    df = pd.DataFrame(protocols_output)
    colr1, colr2 = st.columns([2, 1], gap="small")
    colr1.dataframe(df.style.map(
        highlight_buffer, subset=['NEBuffer r1.1', 'NEBuffer r2.1', 'NEBuffer r3.1', 'rCutSmart'] if
        db == "NEB" else ['Buffer A', 'Buffer B', 'Buffer C', 'Buffer D', 'Buffer E', 'Buffer F', 'Buffer G',
                          'Buffer H', 'Buffer J', 'Buffer K', 'Buffer MultiCore']).apply(
        lambda x: highlight_row(x, db) if x.name in [len(df) - 2, len(df) - 1] else [''] * len(x),
        axis=1), hide_index=True, key=f"{db}_proto")

    if db == "NEB":
        colr2.write("**Protocol**")
        colr2.write("**1.** Prepare mix.")
        df_selected = pd.DataFrame(enzymes)

        unique_temps = df_selected['Incubation Temp (°C)'].unique()
        timesave = True if df_selected['TimeSaver'].eq('✅').all() else False

        if len(unique_temps) == 1:
            colr2.write(f"**2.** Incubate at {unique_temps[0]}°C for {'10min' if timesave is True else '1h'}")
        else:
            colr2.write(f"**2.** Incubation")
            colr2.warning(
                "Selected enzymes have different incubation temperatures. Please choose the appropriate temperature for each enzyme:")
            for temp in unique_temps:
                enzymes_with_temp = df_selected[df_selected['Incubation Temp (°C)'] == temp]['Enzyme'].tolist()
                colr2.write(f"Enzymes with {temp}°C: {', '.join(enzymes_with_temp)}")
        if timesave is False:
            colr2.warning("Some enzymes are not TimeSaver, so we recommend leaving 1 hour or more.")
        else:
            colr2.warning(
                "Even though TimeSaver enzymes are very practical, don't hesitate to add incubation time if necessary.")

        max_heatInactivationTemp = df_selected['heatInactivationTemp'].max()
        max_heatInactivationTime = df_selected['heatInactivationTime'].max()
        if int(max_heatInactivationTemp) > 0:
            colr2.write(
                f"**3.** Inactivate at {max_heatInactivationTemp}°C for {max_heatInactivationTime}min")
            colr2.warning(
                "In my experience, with NEB, if enzymes have different inactivation temperatures and times, it is preferable to use the highest temperature and the longest time.")

        if (df_selected['heatInactivationTemp'] == "0").any():
            colr2.write("**4.** Certain enzymes do not inactivate. A gel purification is required.")
        else:
            colr2.warning("**3. bonus** Gel purification never hurts")
    else:
        colr2.write("**Protocol**")
        colr2.write("**1.** Prepare mix.")
        df_selected = pd.DataFrame(enzymes)

        unique_temps = df_selected['Incubation Temp (°C)'].unique()

        if len(unique_temps) == 1:
            colr2.write(f"**2.** Incubate at {unique_temps[0]}°C for 1h")
        else:
            colr2.write(f"**2.** Incubation")
            colr2.warning(
                "Selected enzymes have different incubation temperatures. Please choose the appropriate temperature for each enzyme:")
            for temp in unique_temps:
                enzymes_with_temp = df_selected[df_selected['Incubation Temp (°C)'] == temp]['Enzyme'].tolist()
                colr2.write(f"- Enzymes with {temp}°C: {', '.join(enzymes_with_temp)}")
        colr2.warning("Don't hesitate to add incubation time if necessary.")

        heatinactivation = True if df_selected['Heat Inactivation'].eq('✅').all() else False
        if heatinactivation is True:
            colr2.write(
                f"**3.** Inactivate at 65°C for 15min")
            colr2.warning(
                "Don't hesitate to add incubation time if necessary. And gel purification never hurts")
        else:
            colr2.write(
                "**3.** Certain enzymes do not inactivate or will damage your DNA. A gel purification is required.")
            non_inactivating_enzymes = df_selected[df_selected['Heat Inactivation'].ne('✅')]

            for idx, enzyme_row in non_inactivating_enzymes.iterrows():
                colr2.write(f"- {enzyme_row['Enzyme']}: Heat Inactivation status - {enzyme_row['Heat Inactivation']}")


def highlight_buffer(val):
    val = str(val)
    color = ''
    if val == "" or val == "nan":
        color = ''
    elif "*" in val or int(val) == 0:
        color = 'background-color: #ff3838;'
    elif 0 < int(val) < 25:
        color = 'background-color: #fb6e16;'
    elif 25 <= int(val) < 50:
        color = 'background-color: #ef9700;'
    elif 50 <= int(val) < 75:
        color = 'background-color: #dcbb00;'
    elif 75 <= int(val) < 100:
        color = 'background-color: #c3db2c;'
    elif int(val) == 100:
        color = 'background-color: #a6f863;'
    return color


def highlight_row(s, db="NEB"):
    vol_tol = 50 if db == "NEB" else 20
    if s.iloc[1] > vol_tol or s.iloc[1] < 0:
        return ['background-color: red'] * len(s)
    else:
        return [''] * len(s)


# Page config
page_config()

NEBtab, PROMEGAtab = st.tabs(["NEB", "Promega"])

with NEBtab:
    enzymes = display_table("NEB")
    select_enzymes = select_enzyme(enzymes, "NEB")
    if select_enzymes:
        digestion_protocols(select_enzymes, "NEB")

with PROMEGAtab:
    enzymes = display_table("PROMEGA")
    select_enzymes = select_enzyme(enzymes, "PROMEGA")
    if select_enzymes:
        digestion_protocols(select_enzymes, "PROMEGA")
