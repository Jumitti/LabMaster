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
import math

import pandas as pd


local_test = platform.processor()

# Settings for Streamlit page
st.set_page_config(page_title="LabMaster", page_icon="🔬", initial_sidebar_state="expanded", layout="wide")

# Main page
st.logo("img/labmaster_logo.png")
st.sidebar.image("img/labmaster_logo.png")
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
st.sidebar.divider()
st.sidebar.page_link("labmaster.py", label="Home")
st.sidebar.page_link("pages/westernblot.py", label="Western Blot")
st.sidebar.page_link("pages/transfection.py", label="Transfection")
st.sidebar.page_link("pages/RT-qPCR.py", label="RT-qPCR")

# MIT licence
st.sidebar.divider()
st.sidebar.link_button('Under MIT licence', 'https://github.com/Jumitti/BlotMaster/blob/main/LICENSE')

# Main
st.markdown("<h3 style='text-align: center; color: black;'>Welcome to LabMaster</h1>", unsafe_allow_html=True)
st.markdown('')

st.image('img/labmaster_banner.png')
st.markdown('')

st.markdown('**Overview**')
st.markdown(
        '<div style="text-align: justify;">LabMaster is designed to streamline and simplify the calculations essential for biological research. Our software currently offers tools for transfection calculations and sample preparation for Western blot analysis. Here is a brief overview of the features:</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/transfection.py", label='**Transfection Calculations**',)
st.markdown(
        '<div style="text-align: justify;">Our transfection module assists researchers in accurately determining the necessary quantities of DNA, reagents, and other components required for efficient and effective transfection. Whether you are working with plasmid DNA, siRNA, or other nucleic acids, LabMaster ensures precise calculations to optimize your experimental outcomes.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/westernblot.py", label="**Western Blot Sample Preparation**")
st.markdown(
        '<div style="text-align: justify;">The Western blot module facilitates the preparation of samples for protein analysis. This tool helps you calculate the exact volumes and concentrations needed for your samples, ensuring reproducibility and consistency in your experiments. From lysate preparation to loading buffer calculations, LabMaster supports every step of the process.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.markdown('**Community Project**')
st.markdown(
        '<div style="text-align: justify;">LabMaster is a community-driven project. We welcome contributions and suggestions from all users. Whether you want to propose new features, improve existing tools, or contribute code, your input is valued and appreciated. Together, we can enhance LabMaster to better serve the needs of the biological research community.</div>',
        unsafe_allow_html=True)
st.markdown(
        '<div style="text-align: justify;">We are committed to expanding LabMaster with additional features to cover more areas of biological research, making your lab work more efficient and less error-prone.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.markdown('**Disclaimer**')
st.markdown(
        '<div style="text-align: justify;">LabMaster is intended to serve as an aid for researchers, providing assistance with common laboratory calculations. While we strive for accuracy, we cannot guarantee that all calculations will be free from errors. Users are advised to double-check critical calculations and use their own scientific judgment in conjunction with this tool.</div>',
        unsafe_allow_html=True)
st.markdown("")
