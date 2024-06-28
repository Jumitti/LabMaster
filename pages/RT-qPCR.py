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

import pandas as pd
import streamlit as st

from utils.page_config import page_config

# Page config
page_config()

# Sample table
colm1, colm2, colm3, colm4 = st.columns([1, 0.5, 1, 1], gap="small")
colm1.write("**Input user**")

initial_data = [
    {"Sample name": "WT_1", "Conc. 1 (ng/µL)": 253.5, "Conc. 2 (ng/µL)": 256},
    {"Sample name": "KO_1", "Conc. 1 (ng/µL)": 300, "Conc. 2 (ng/µL)": 302.2}
]

if "samples_table_rna" not in st.session_state:
    st.session_state["samples_table_rna"] = pd.DataFrame(initial_data)

samples_table = colm1.data_editor(st.session_state["samples_table_rna"],
                                       num_rows="dynamic", hide_index=True)

if colm1.button("Save"):
    st.session_state["samples_table_rna"] = samples_table

if "Conc. 1 (ng/µL)" in samples_table.columns and "Conc. 2 (ng/µL)" in samples_table.columns:
    samples_table["Mean (ng/µL)"] = samples_table[["Conc. 1 (ng/µL)", "Conc. 2 (ng/µL)"]].mean(axis=1)
elif "Conc. 1 (ng/µL)" in samples_table.columns:
    samples_table["Mean (ng/µL)"] = samples_table["Conc. 1 (ng/µL)"].mean(axis=1)
elif "Conc. 2 (ng/µL)" in samples_table.columns:
    samples_table["Mean (ng/µL)"] = samples_table["Conc. 2 (ng/µL)"].mean(axis=1)
else:
    samples_table["Mean (ng/µL)"] = 0

colm2.write("**Summary and samples average**")
mean_samples_table = samples_table[["Sample name", "Mean (ng/µL)"]]
colm2.dataframe(mean_samples_table, hide_index=True)

colm3.write("**Samples preparation**")
elution_volume = colm3.number_input("Elution volume (µL)", value=30.00, min_value=0.01, step=0.01)
experimental_rna = colm3.number_input("Experimental RNA (ng)", value=2000.00, min_value=0.01, step=0.01)
primer = colm3.number_input("Primer (µL)", value=1.00, min_value=0.01, step=0.01)
final_volume = colm3.number_input("Final volume (µL)", value=10, min_value=5, step=5)

colm4.write("**After RT settings**")
final_concentration = colm4.number_input("Final concentration (ng/µL)", value=20.00, min_value=0.01, step=0.01)

st.divider()

def sample_vol(row):
    return experimental_rna / row["Mean (ng/µL)"]

def h2o_after_rt(row):
    if row["H2O (µL)"] >= 0:
        return (experimental_rna / final_concentration) - (2 * final_volume)
    else:
        return (experimental_rna / final_concentration) - (2 * final_volume) + row["H2O (µL)"]


# Output Table for Mix
required_columns = ["Sample name", "Mean (ng/µL)"]
if not set(required_columns).issubset(mean_samples_table.columns):
    st.error("The columns 'Sample name' and 'Mean (ng/µL)' are required.")
else:
    missing_data = mean_samples_table[
        mean_samples_table["Sample name"].isnull() | mean_samples_table["Mean (ng/µL)"].isnull()]
    if not missing_data.empty or (mean_samples_table["Mean (ng/µL)"] < 0).any():
        st.error("Some lines do not have concentrations or are negative.")
    else:
        adjusted_samples_table = mean_samples_table.copy()
        adjusted_samples_table["Sample volume (µL)"] = adjusted_samples_table.apply(sample_vol, axis=1)
        adjusted_samples_table["Primer (µL)"] = primer
        adjusted_samples_table["H2O (µL)"] = final_volume - adjusted_samples_table["Sample volume (µL)"] - primer
        adjusted_samples_table["RT mix (µL)"] = final_volume
        adjusted_samples_table["H2O after RT (µL)"] = adjusted_samples_table.apply(h2o_after_rt, axis=1)

        st.dataframe(adjusted_samples_table, hide_index=True)

        RT_mix = [
            {"Component": "GoScript 5X Reaction Buffer", "Volume (µL)": 2 * final_volume/5},
            {"Component": "MgCl2 (25 mM)", "Volume (µL)": final_volume * 1.6 / 10},
            {"Component": "dNTP (10 mM)", "Volume (µL)": final_volume / 10},
            {"Component": "GoScript Reverse Transcriptase", "Volume (µL)": final_volume / 10},
            {"Component": "H2O", "Volume (µL)": final_volume * 0.24}
        ]
        RT_mix_table = pd.DataFrame(RT_mix)
        st.dataframe(RT_mix_table, hide_index=True)
