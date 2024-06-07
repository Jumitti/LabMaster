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

import streamlit as st
import platform

import pandas as pd


local_test = platform.processor()

# Settings for Streamlit page
st.set_page_config(page_title="BlotMaster", page_icon="🔬", initial_sidebar_state="expanded", layout="wide")

# Main page
st.sidebar.title('👩🏼‍🔬 BlotMaster')
st.sidebar.write("Created by Minniti Julien")

# Button sidebar
colsb1, colsb2, colsb3 = st.sidebar.columns(3, gap="small")
colsb1.link_button("Help ⁉",
                  '')
colsb2.link_button('GitHub', 'https://github.com/Jumitti/BlotMaster')
if local_test == "":
    colsb3.link_button('Download app 📥', 'https://github.com/Jumitti/BlotMaster/releases')
else:
    colsb3.link_button('Web app 🌐', 'https://blotmaster.streamlit.app/')

# MIT licence
st.sidebar.divider()
st.sidebar.link_button('Under MIT licence', 'https://github.com/Jumitti/BlotMaster/blob/main/LICENSE')

# Main
# Sample table
colm1, colm2, colm3 = st.columns([1, 1, 1], gap="small")
colm1.write("**Samples table**")
df = pd.DataFrame(
    [{"Sample name": "WT_1", "µg/µL": 1.5}, {"Sample name": "KO_1", "µg/µL": 2}]
)
df = colm1.data_editor(df, key="vector", num_rows="dynamic")

dilution_factor = colm2.number_input("Dilution factor of samples:", min_value=1.0, step=0.1, help='e.g. if you measured the proteins in 2 µL of your sample with a Bradford or other then the dilution factor of your sample is 2')
sample_volume = colm2.number_input("Total sample volumes (µL):", min_value=0.00, step=0.01)

proteins_per_well = colm3.number_input("Amount of proteins/well (µg):", min_value=0.00, step=0.01)
volume_per_well = colm3.number_input("Total volume/well (µL):", min_value=0.00, step=0.01, help='We take into account the sample, the lysis buffer and the load buffer')
concentration_samples = proteins_per_well/volume_per_well if volume_per_well > 0 else proteins_per_well
colm3.write(f"Concentration of samples (µg/µL): {concentration_samples}")
st.divider()


