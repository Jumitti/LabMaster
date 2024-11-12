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
import math
from datetime import datetime

import pandas as pd
import streamlit as st

from utils.page_config import page_config


# For Master Mix
def adjust_concentration(row):
    return row["Mean (µg/µL)"] / dilution_factor


def LyB_to_add(row):
    return (((row["Adjusted conc. (µg/µL)"] * sample_volume) / concentration_samples) - sample_volume) - (
            (row["Adjusted conc. (µg/µL)"] * sample_volume) / concentration_samples) / LoB_concentration


def loading_buffer(row):
    if row["Lysis buffer to add (µL)"] <= 0:
        return sample_volume / (LoB_concentration - 1)
    else:
        return ((row["Adjusted conc. (µg/µL)"] * sample_volume) / concentration_samples) / LoB_concentration


def total_volume(row):
    if row["Lysis buffer to add (µL)"] >= 0:
        return row["Lysis buffer to add (µL)"] + sample_volume + row["Loading buffer (µL)"]
    else:
        return sample_volume + row["Loading buffer (µL)"]


def vol_to_charge(row):
    if row["Adjusted conc. (µg/µL)"] * sample_volume < proteins_per_well:
        return "Proteins level too low"
    elif row["Lysis buffer to add (µL)"] <= 0:
        return row["Total volume (µL)"]
    else:
        return float(volume_per_well)


def nb_sample(row):
    if row["Adjusted conc. (µg/µL)"] * sample_volume < proteins_per_well:
        return 0
    elif row["Lysis buffer to add (µL)"] <= 0:
        return 1
    else:
        return math.floor((row["Total volume (µL)"]) / volume_per_well)


# For one sample to load
def sample_vol_OS(row):
    return proteins_per_well / row["Adjusted conc. (µg/µL)"]


def LyB_to_add_OS(row):
    return volume_per_well - row["Sample volume (µL)"] - volume_per_well / LoB_concentration


def loading_buffer_OS(row):
    if row["Lysis buffer to add (µL)"] < 0:
        return row["Sample volume (µL)"] / (LoB_concentration - 1)
    else:
        return volume_per_well / LoB_concentration


def total_volume_OS(row):
    if row["Lysis buffer to add (µL)"] >= 0:
        return row["Lysis buffer to add (µL)"] + row["Sample volume (µL)"] + row["Loading buffer (µL)"]
    else:
        return row["Sample volume (µL)"] + row["Loading buffer (µL)"]


def vol_to_charge_OS(row):
    if row["Adjusted conc. (µg/µL)"] * sample_volume < proteins_per_well:
        return "Proteins level too low"
    elif row["Lysis buffer to add (µL)"] <= 0:
        return row["Total volume (µL)"]
    else:
        return float(volume_per_well)


def nb_sample_OS(row):
    if row["Adjusted conc. (µg/µL)"] * sample_volume < proteins_per_well:
        return 0
    elif row["Lysis buffer to add (µL)"] <= 0:
        return 1
    else:
        return math.floor(sample_volume / row["Sample volume (µL)"])


# Export to Excel
def to_excel_with_style(df1, df2, df3):
    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')

    df1.to_excel(writer, index=False, sheet_name='One sample')
    workbook = writer.book
    worksheet1 = writer.sheets['One sample']
    red_format = workbook.add_format({'bg_color': '#FF0000'})
    orange_format = workbook.add_format({'bg_color': '#FFA500'})
    apply_formatting(worksheet1, df1, red_format, orange_format)
    df2.to_excel(writer, index=False, sheet_name='Master mix')
    worksheet2 = writer.sheets['Master mix']
    apply_formatting(worksheet2, df2, red_format, orange_format)
    df3.to_excel(writer, index=False, sheet_name='Settings')

    writer.close()
    processed_data = output.getvalue()
    return processed_data


def apply_formatting(worksheet, df, red_format, orange_format):
    for row_num, value in enumerate(df['Nb of samples'], start=1):
        cell_format = None
        if value <= 1:
            cell_format = red_format
        elif value == 2:
            cell_format = orange_format
        if cell_format:
            worksheet.set_row(row_num, None, cell_format)


# Page config
page_config()

# Sample table
colm1, colm2, colm3, colm4 = st.columns([1, 0.5, 1, 1], gap="small")
colm1.write("**Input user**")

initial_data = [
    {"Sample name": "WT_1", "Conc. 1 (µg/µL)": 1.5, "Conc. 2 (µg/µL)": 2.0},
    {"Sample name": "KO_1", "Conc. 1 (µg/µL)": 2.0, "Conc. 2 (µg/µL)": 3.0}
]
if "samples_table" not in st.session_state:
    st.session_state["samples_table"] = pd.DataFrame(initial_data)
elif "Mean (µg/µL)" in st.session_state["samples_table"].columns:
    st.session_state["samples_table"] = st.session_state["samples_table"].drop(columns=["Mean (µg/µL)"])

samples_table = colm1.data_editor(
    st.session_state["samples_table"],
    num_rows="dynamic",
    hide_index=True
)

if colm1.button("Save"):
    st.session_state["samples_table"] = samples_table

if "Conc. 1 (µg/µL)" in samples_table.columns and "Conc. 2 (µg/µL)" in samples_table.columns:
    samples_table["Mean (µg/µL)"] = samples_table[["Conc. 1 (µg/µL)", "Conc. 2 (µg/µL)"]].mean(axis=1)
elif "Conc. 1 (µg/µL)" in samples_table.columns:
    samples_table["Mean (µg/µL)"] = samples_table["Conc. 1 (µg/µL)"]
elif "Conc. 2 (µg/µL)" in samples_table.columns:
    samples_table["Mean (µg/µL)"] = samples_table["Conc. 2 (µg/µL)"]
else:
    samples_table["Mean (µg/µL)"] = 0

colm2.write("**Summary and samples average**")
mean_samples_table = samples_table[["Sample name", "Mean (µg/µL)"]]
colm2.dataframe(mean_samples_table, hide_index=True)

# Settings
dilution_factor = colm3.number_input("**Dilution factor of samples:**", min_value=0.01, step=0.01,
                                     value=1.00 if "dilution_factor_save" not in st.session_state else st.session_state['dilution_factor_save'],
                                     help='e.g. if you measured the proteins in 2 µL of your sample with a Bradford or '
                                          'other then the dilution factor of your sample is 2')
st.session_state["dilution_factor_save"] = dilution_factor

sample_volume = colm3.number_input("**Total sample volumes (µL):**", min_value=0.01, step=0.01,
                                   value=30.00 if "sample_volume_save" not in st.session_state else st.session_state['sample_volume_save'])
st.session_state["sample_volume_save"] = sample_volume

proteins_per_well = colm4.number_input("**Amount of proteins/well (µg):**", min_value=0.01, step=0.01,
                                       value=25.00 if "proteins_per_well_save" not in st.session_state else st.session_state['proteins_per_well_save'])
st.session_state["proteins_per_well_save"] = proteins_per_well

volume_per_well = colm4.number_input("**Total volume/well (µL):**", min_value=0.01, step=0.01,
                                     value=30.00 if "volume_per_well_save" not in st.session_state else st.session_state['volume_per_well_save'],
                                     help='We take into account the sample, the lysis buffer and the load buffer')
st.session_state["volume_per_well_save"] = volume_per_well

concentration_samples = proteins_per_well / volume_per_well
colm4.write(f"Concentration of samples: {concentration_samples} µg/µL")

LoB_concentration = colm4.number_input("**Loading buffer conc. (X):**", min_value=2.00, step=0.01,
                                       value=5.00 if "LoB_concentration_save" not in st.session_state else st.session_state['LoB_concentration_save'],
                                       help="e.g. the loading buffer is often in concentration 'X'. If it's 5X then put 5")
st.session_state["LoB_concentration_save"] = LoB_concentration

settings = [{"Dilution factor of samples": dilution_factor, "Total sample volumes (µL)": sample_volume,
             "Amount of proteins/well (µg)": proteins_per_well, "Total volume/well (µL)": volume_per_well,
             "Concentration of samples (µg/µL)": concentration_samples, "Loading buffer conc. (X):": LoB_concentration}]
settings = pd.DataFrame(settings)
st.divider()

# Output Table for Mix
required_columns = ["Sample name", "Mean (µg/µL)"]
if not set(required_columns).issubset(mean_samples_table.columns):
    st.error("The columns 'Sample name' and 'Mean (µg/µL)' are required.")
else:
    missing_data = mean_samples_table[
        mean_samples_table["Sample name"].isnull() | mean_samples_table["Mean (µg/µL)"].isnull()]
    if not missing_data.empty or (mean_samples_table["Mean (µg/µL)"] <= 0).any():
        st.error("Some lines do not have concentrations or are negative.")
    else:
        adjusted_samples_table = mean_samples_table.copy()
        adjusted_samples_table["Adjusted conc. (µg/µL)"] = adjusted_samples_table.apply(adjust_concentration, axis=1)
        adjusted_samples_table["Lysis buffer to add (µL)"] = adjusted_samples_table.apply(LyB_to_add, axis=1)
        adjusted_samples_table["Loading buffer (µL)"] = adjusted_samples_table.apply(loading_buffer, axis=1)
        adjusted_samples_table["Total volume (µL)"] = adjusted_samples_table.apply(total_volume, axis=1)
        adjusted_samples_table["Volume to charge (µL)"] = adjusted_samples_table.apply(vol_to_charge, axis=1)
        adjusted_samples_table["Nb of samples"] = adjusted_samples_table.apply(nb_sample, axis=1)

        styled_table = adjusted_samples_table.style.apply(lambda row: [
            'background-color: red' if row["Nb of samples"] <= 1 else 'background-color: orange' if row["Nb of samples"] == 2
            else '' for _ in row], axis=1)

        st.markdown("**For a Master Mix**", help="Here, the sample is prepared for several samples at the same time, "
                                                 "depending on the volume of char in the Western Blot (it's a Master Mix).")
        st.dataframe(styled_table, hide_index=True)

        one_sample = mean_samples_table.copy()
        one_sample["Adjusted conc. (µg/µL)"] = adjusted_samples_table["Adjusted conc. (µg/µL)"]
        one_sample["Sample volume (µL)"] = one_sample.apply(sample_vol_OS, axis=1)
        one_sample["Lysis buffer to add (µL)"] = one_sample.apply(LyB_to_add_OS, axis=1)
        one_sample["Loading buffer (µL)"] = one_sample.apply(loading_buffer_OS, axis=1)
        one_sample["Total volume (µL)"] = one_sample.apply(total_volume_OS, axis=1)
        one_sample["Volume to charge (µL)"] = one_sample.apply(vol_to_charge_OS, axis=1)
        one_sample["Nb of samples"] = one_sample.apply(nb_sample_OS, axis=1)

        styled_table_OS = one_sample.style.apply(lambda row: [
            'background-color: red' if row["Nb of samples"] <= 1 else 'background-color: orange' if row["Nb of samples"] == 2
            else '' for _ in row], axis=1)

        st.markdown("**For one sample**",
                    help="Here, the sample is prepared for several samples at the same time, depending on the volume of "
                         "char in the Western Blot (it's a Master Mix).")
        st.dataframe(styled_table_OS, hide_index=True)

        df_xlsx = to_excel_with_style(one_sample, adjusted_samples_table, settings)
        current_date_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        st.download_button(label='💾 Download tables (.xlsx)',
                           data=df_xlsx,
                           file_name=f'WB_samples_{current_date_time}.xlsx',
                           mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
