# MIT License
#
# Copyright (c) 2024 Minniti Julien
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import io
import platform
from datetime import datetime

import pandas as pd
import streamlit as st

local_test = platform.processor()

culture_vessel_options_lipo = {"96-well": 5, "24-well": 25, "12-well": 50, "6-well/35-mm": 125,
                               "60-mm/flask 25 cm2": 250, "10-cm/flask 75 cm2": 500, "15-cm/flask 175cm2": 750}
culture_vessel_options_jetprime = {"96-well": 5, "24-well": 50, "12-well": 75, "6-well/35-mm": 200,
                                   "60-mm/flask 25 cm2": 200, "10-cm/flask 75 cm2": 500, "15-cm/flask 175cm2": 1000}
calcul = []
calcul_j = []

local_test = platform.processor()

# Settings for Streamlit page
st.set_page_config(page_title="LabMaster", page_icon="🔬", initial_sidebar_state="expanded", layout="wide")

# Main page
st.sidebar.title('👩🏼‍🔬 LabMaster')
st.sidebar.write("Created by Minniti Julien")

# Button sidebar
colsb1, colsb2, colsb3 = st.sidebar.columns(3, gap="small")
colsb1.link_button("Help ⁉", '')
colsb2.link_button('GitHub', 'https://github.com/Jumitti/BlotMaster')
if local_test == "":
    colsb3.link_button('Download app 📥', 'https://github.com/Jumitti/BlotMaster/releases')
else:
    colsb3.link_button('Web app 🌐', 'https://blotmaster.streamlit.app/')

# Table
st.sidebar.page_link("labmaster.py", label="Main")
st.sidebar.page_link("pages/westernblot.py", label="Western Blot")
st.sidebar.page_link("pages/transfection.py", label="Transfection")

# Settings for all transfection
st.sidebar.divider()
df = pd.DataFrame(
    [{"Plasmid/Vector/RNA": "pCDNA3.1", "µg/µL": 1.5}, {"Plasmid/Vector/RNA": "pSBBi-GP", "µg/µL": 2}]
)
st.sidebar.write("**List of Plasmid/Vector/RNA**")
df = st.sidebar.data_editor(df, key="vector", num_rows="dynamic")
vector_table = df.to_dict(orient='records')

st.sidebar.write("**Number of condition(s)**")
num_columns = st.sidebar.number_input('Number of condition(s)', min_value=1, max_value=None, value=1, step=1, help=None,
                                      key="mix", label_visibility="collapsed")

st.sidebar.divider()
st.sidebar.write("**Settings for all conditions**")
vector_for_all_mix = st.sidebar.toggle("Vector(s) per condition")
if vector_for_all_mix:
    vector_per_mix = st.sidebar.number_input('Vector per condition', min_value=1, max_value=len(df), value=1, step=1,
                                             help=None,
                                             key=f"vector_per_mix")

amount_dna_for_all = st.sidebar.toggle("Amount of DNA per well",
                                       help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
if amount_dna_for_all:
    amount_of_dna = st.sidebar.number_input('Amount of DNA per well (µg)', min_value=0.00, max_value=30.00, value=1.00,
                                            step=0.01, key=f'amount_dna')

selected_vector = False
if vector_for_all_mix:
    selected_vector = st.sidebar.toggle('Selection of vector')
    if selected_vector:
        for j in range(1, vector_per_mix + 1):
            vector_selected = st.sidebar.selectbox(f'Vector {j}', df["Plasmid/Vector/RNA"], key=f'vector_selected{j}')
            if amount_dna_for_all:
                sum_amount_dna_selected_same_i = sum(item[2] for item in calcul_j)
                max_slider = amount_of_dna - sum_amount_dna_selected_same_i
                if max_slider > 0:
                    amount_dna_selected = st.sidebar.slider(
                        f'Amount of {vector_selected if vector_selected is not None or vector_selected != "" else f"Vector {j}"} (µg)',
                        0.00, max_slider, value=max_slider, step=0.01, key=f'amount_vector_selected{j}')
                else:
                    st.sidebar.warning(
                        f"Amount of {vector_selected if vector_selected is not None or vector_selected != '' else f'Vector {j}'}"
                        f"not available. The previous plasmid(s)/vector(s)/RNA already uses all the DNA required for transfection.")
                    amount_dna_selected = 0
            else:
                st.sidebar.warning(
                    "Amount of vector available only if you select the option to define the amount of DNA desired")
            calcul_j.append([j, vector_selected, amount_dna_selected if amount_dna_for_all else ""])

transfection_type_for_all = st.sidebar.toggle("Transfection type",
                                              help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
if transfection_type_for_all:
    transfection_type = st.sidebar.radio("Transfection type", ["Lipofectamine (2000/3000)", 'jetPRIME'],
                                         key=f"transfection_type", horizontal=True)
    if transfection_type == "Lipofectamine (2000/3000)":
        transfection_type_index = 0
    else:
        transfection_type_index = 1
else:
    transfection_type_index = 0

dna_ratio_for_all = st.sidebar.toggle("DNA ratio",
                                      help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
if dna_ratio_for_all:
    st.sidebar.write("**Ratio DNA(µg)/Reagent(µL)**")
    col1, col2 = st.sidebar.columns(2, gap="small")
    dna_ratio = col1.number_input("DNA", min_value=0.1, max_value=None, value=1.0, step=0.1, help=None,
                                  key=f"dna_ratio")
    reagent_ratio = col2.number_input("Reagent", min_value=0.1, max_value=None, value=2.0, step=0.1,
                                      help=None, key=f"reagent_ratio")

culture_vessel_for_all = st.sidebar.toggle('Culture Vessel')
if culture_vessel_for_all:
    culture_vessel = st.sidebar.selectbox("Culture Vessel", list(culture_vessel_options_lipo.keys()),
                                          key=f"culture_vessel")

numbers_wells_for_all = st.sidebar.toggle("Number of well(s)")
if numbers_wells_for_all:
    number_wells = st.sidebar.number_input("Number of well(s)", min_value=1, max_value=None, value=1, step=1, help=None,
                                           key=f"number_well")

# MIT licence
st.sidebar.divider()
st.sidebar.link_button('Under MIT licence', 'https://github.com/Jumitti/BlotMaster/blob/main/LICENSE')

columns = st.columns(num_columns)

try:
    for i in range(1, num_columns + 1):
        with columns[i % num_columns - 1]:
            st.write('**Condition name**')
            name_condition = st.text_input(f"Condition name", value=f"Condition  {i}", max_chars=20,
                                           key=f"condition_name{i}",
                                           label_visibility="collapsed")

            st.write(f'**Number of vector(s) in {name_condition}**')
            vector_per_mix = st.number_input('Vector(s) per condition', min_value=1, max_value=len(df),
                                             value=1 if vector_for_all_mix is False else vector_per_mix, step=1,
                                             help=None, key=f"vector_per_mix{i}", label_visibility="collapsed")

            st.markdown(f'**Amount of DNA per well (µg)**',
                        help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
            amount_of_dna = st.number_input('Amount of DNA per well (µg)', min_value=0.00, max_value=30.00,
                                            value=1.00 if amount_dna_for_all is False else amount_of_dna, step=0.01,
                                            help=None, key=f'amount_dna{i}', label_visibility="collapsed")

            st.markdown(f'**Transfection type**',
                        help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
            transfection_type = st.radio("Transfection type", ["Lipofectamine (2000/3000)", 'jetPRIME'],
                                         index=transfection_type_index, key=f"transfection_type{i}", horizontal=True,
                                         label_visibility="collapsed")

            st.markdown("**Ratio DNA(µg)/Reagent(µL)**",
                        help=f"Please, see [HELP](https://www.notion.so/jumitti/Transfection-Ratio-Xpert-c25f10537e454529a4397925fe9e68c5)")
            col1, col2 = st.columns(2)
            dna_ratio = col1.number_input("DNA", min_value=0.1, max_value=None,
                                          value=1.0 if dna_ratio_for_all is False else dna_ratio, step=0.1, help=None,
                                          key=f"dna_ratio{i}")
            reagent_ratio = col2.number_input("Reagent", min_value=0.1, max_value=None,
                                              value=2.0 if dna_ratio_for_all is False else reagent_ratio, step=0.1,
                                              help=None, key=f"reagent_ratio{i}")
            ratio = dna_ratio / reagent_ratio

            col1b, col2b = st.columns(2)
            col1b.write("**Culture Vessel**")
            culture_vessel = col1b.selectbox("Culture Vessel", list(culture_vessel_options_lipo.keys()),
                                             index=list(culture_vessel_options_lipo.keys()).index(st.session_state.get("culture_vessel")) if culture_vessel_for_all else 0,
                                             key=f"culture_vessel{i}", label_visibility="collapsed")

            col2b.write("**Number of well(s)**")
            number_wells = col2b.number_input("Number of well(s)", min_value=1, max_value=None,
                                              value=1 if numbers_wells_for_all is False else number_wells, step=1,
                                              help=None, key=f"number_well{i}", label_visibility="collapsed")
            culture_vessel_value = culture_vessel_options_lipo[
                culture_vessel] if transfection_type == "Lipofectamine (2000/3000)" \
                else culture_vessel_options_jetprime[culture_vessel]

            for j in range(1, vector_per_mix + 1):
                st.write(f'**Selection of vector {j}**')
                vector_selected = st.selectbox(f'Vector {j}',
                                               df["Plasmid/Vector/RNA"],
                                               index=df["Plasmid/Vector/RNA"].to_list().index(st.session_state.get(
                                                   f"vector_selected{j}")) if selected_vector else 0,
                                               key=f'vector_selected{i}-{j}', label_visibility="collapsed")

                sum_amount_dna_selected_same_i = sum(item[8] for item in calcul if item[0] == i)
                max_slider = amount_of_dna - sum_amount_dna_selected_same_i

                if max_slider > 0:
                    st.write(
                        f'**Amount of {vector_selected if vector_selected is not None or vector_selected != "" else f"Vector {j}"} (µg)**')
                    amount_dna_selected = st.slider(f'DNA', 0.00, max_slider,
                                                    value=max_slider if selected_vector is False or amount_dna_for_all is False
                                                    else st.session_state[f"amount_vector_selected{j}"], step=0.01,
                                                    key=f'amount_vector_selected{i}-{j}', label_visibility="collapsed")
                else:
                    st.warning(
                        f"Amount of {vector_selected if vector_selected is not None or vector_selected != '' else f'Vector {j}'}"
                        f"not available. The previous plasmid(s)/vector(s)/RNA already uses all the DNA required for transfection.")
                    amount_dna_selected = 0
                calcul.append(
                    [i, name_condition, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells,
                     vector_selected, amount_dna_selected])

    results_dict = {}
    for i, name_condition, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells, vector_selected, amount_dna_selected in calcul:
        if vector_selected not in results_dict:
            results_dict[vector_selected] = {
                "Plasmid/Vector/RNA": vector_selected,
                f"{name_condition} (µL)": (amount_dna_selected / next(item['µg/µL'] for item in vector_table
                                                                      if item[
                                                                          "Plasmid/Vector/RNA"] == vector_selected)) * number_wells
            }
        else:
            results_dict[vector_selected][f"{name_condition} (µL)"] = (amount_dna_selected / next(
                item['µg/µL'] for item in vector_table
                if item["Plasmid/Vector/RNA"] == vector_selected)) * number_wells

        if transfection_type not in results_dict:
            results_dict[transfection_type] = {
                "Plasmid/Vector/RNA": transfection_type,
                f"{name_condition} (µL)": (amount_of_dna / ratio) * number_wells}
        else:
            results_dict[transfection_type][f"{name_condition} (µL)"] = (amount_of_dna / ratio) * number_wells

        if transfection_type == "Lipofectamine (2000/3000)":
            if "Culture Vessel Lipo 1" not in results_dict:
                results_dict["Culture Vessel Lipo 1"] = {
                    "Plasmid/Vector/RNA": "OptiMEM w/ Plasmid/Vector/RNA",
                    f"{name_condition} (µL)": culture_vessel_value * number_wells}
                results_dict["Culture Vessel Lipo 2"] = {
                    "Plasmid/Vector/RNA": "OptiMEM w/ Lipofectamine",
                    f"{name_condition} (µL)": culture_vessel_value * number_wells}
            else:
                results_dict["Culture Vessel Lipo 1"][f"{name_condition} (µL)"] = culture_vessel_value * number_wells
                results_dict["Culture Vessel Lipo 2"][f"{name_condition} (µL)"] = culture_vessel_value * number_wells

        elif transfection_type == "jetPRIME":
            if "Culture Vessel JP" not in results_dict:
                results_dict["Culture Vessel JP"] = {"Plasmid/Vector/RNA": "jetPRIME Buffer",
                                                     f"{name_condition} (µL)": culture_vessel_value * number_wells}
            else:
                results_dict["Culture Vessel JP"][f"{name_condition} (µL)"] = culture_vessel_value * number_wells

        if "Volume per well" not in results_dict:
            results_dict["Volume per well"] = {"Plasmid/Vector/RNA": "Volume per well",
                                               f"{name_condition} (µL)": culture_vessel_value if transfection_type == "jetPRIME"
                                               else culture_vessel_value * 2}
        else:
            results_dict["Volume per well"][
                f"{name_condition} (µL)"] = culture_vessel_value if transfection_type == "jetPRIME" \
                else culture_vessel_value * 2

    for i, name_condition, amount_of_dna, transfection_type, ratio, culture_vessel_value, number_wells, vector_selected, amount_dna_selected in calcul:
        if transfection_type == "Lipofectamine (2000/3000)":
            if "Culture Vessel Lipo 1" in results_dict:
                culture_entry = results_dict.pop("Culture Vessel Lipo 1")
                results_dict["Culture Vessel Lipo 1"] = culture_entry

            if "Culture Vessel Lipo 2" in results_dict:
                culture_lipo_entry = results_dict.pop("Culture Vessel Lipo 2")
                results_dict["Culture Vessel Lipo 2"] = culture_lipo_entry

            if transfection_type in results_dict:
                reagent_entry = results_dict.pop(transfection_type)
                results_dict[transfection_type] = reagent_entry

        else:
            if "Culture Vessel JP" in results_dict:
                culture_entry = results_dict.pop("Culture Vessel JP")
                results_dict["Culture Vessel JP"] = culture_entry

            if transfection_type in results_dict:
                reagent_entry = results_dict.pop(transfection_type)
                results_dict[transfection_type] = reagent_entry

        if "Volume per well" in results_dict:
            vol_entry = results_dict.pop("Volume per well")
            results_dict["Volume per well"] = vol_entry

    df_results = pd.DataFrame(list(results_dict.values()))

    st.dataframe(df_results, hide_index=True)
    csv_file = df_results.to_csv(index=False)
    csv_file2 = df.to_csv(index=False)
    excel_file = io.BytesIO()
    with pd.ExcelWriter(excel_file, engine='xlsxwriter') as writer:
        df_results.to_excel(writer, index=False, sheet_name='Transfection table')
        df.to_excel(writer, index=False, sheet_name='DNA library')

    excel_file.seek(0)
    current_date_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.download_button("💾 Download table (.xlsx)", excel_file,
                       file_name=f'TRXpert_{current_date_time}.xlsx',
                       mime="application/vnd.ms-excel", key='download-excel')

except Exception as e:
    st.error(
        f"Something wrong happened... please contact minnitijulien06@gmail.com of "
        f"[Issues](https://github.com/Jumitti/TransfectionRatioXpert) on GitHub.\n\n {e}")
