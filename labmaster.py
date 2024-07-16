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

from utils.page_config import page_config

# Page config
page_config(logo=True)

# Main
st.markdown("<h3 style='text-align: center; color: black;'>Welcome to LabMaster</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: black;'>👩‍🔬👨‍🔬👩🏻‍🔬👨🏻‍🔬👩🏼‍🔬👨🏼‍🔬👩🏽‍🔬👨🏽‍🔬👩🏾‍🔬👨🏾‍🔬👩🏿‍🔬👨🏿‍🔬</h1>", unsafe_allow_html=True)
st.markdown('')

st.markdown('**Overview**')
st.markdown(
        '<div style="text-align: justify;">LabMaster is designed to streamline and simplify the calculations essential for biological research. Our software currently offers tools for transfection calculations, sample preparation for Western blot analysis, reverse transcription calculations for RNA work, Venn diagram generation, and now includes a Restriction Enzymes tool. Here is a brief overview of the features:</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/restriction_enzyme.py", label='**Restriction Enzymes**',)
st.markdown(
        '<div style="text-align: justify;">With the Restriction Enzymes tool, you can view the specificity of enzymes from NEB and Promega (please note, we are not affiliated with these companies) and generate an enzymatic digestion protocol for molecular biology. This tool helps you determine the best conditions for your restriction enzyme digests, ensuring accurate and reproducible results in your cloning and molecular biology experiments.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/RT-qPCR.py", label='**Reverse Transcription Calculations**',)
st.markdown(
        '<div style="text-align: justify;">Our reverse transcription module provides essential calculations for converting RNA into cDNA. This tool ensures accurate reagent volumes and reaction conditions for efficient reverse transcription, a crucial step in many molecular biology workflows involving RNA.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/transfection.py", label='**Transfection Calculations**',)
st.markdown(
        '<div style="text-align: justify;">Our transfection module assists researchers in accurately determining the necessary quantities of DNA, reagents, and other components required for efficient and effective transfection. Whether you are working with plasmid DNA, siRNA, or other nucleic acids, LabMaster ensures precise calculations to optimize your experimental outcomes.</div>',
        unsafe_allow_html=True)
st.markdown("")

st.page_link("pages/venn_diagram.py", label='**Venn Diagram Generation**',)
st.markdown(
        '<div style="text-align: justify;">LabMaster includes a Venn diagram tool for visualizing the overlap between different data sets. This feature is particularly useful for comparing gene or protein expression profiles, identifying common elements, and presenting your data in an easily interpretable format.</div>',
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

st.markdown('**Contributing**')
st.markdown(
    '''
    <div style="text-align: justify;">
        We welcome contributions from the community! Please read our 
        <a href="https://github.com/Jumitti/LabMaster/blob/main/CONTRIBUTING.md" target="_blank">contributing guidelines</a> 
        to get started. Whether you are fixing bugs, adding new features, or 
        improving documentation, your help is appreciated.
    </div>
    ''',
    unsafe_allow_html=True
)
st.markdown("")

st.markdown('**Contact**')
st.markdown(
    '''
    <div style="text-align: justify;">
        If you have any questions, suggestions, or need support, please reach out to us at 
        <a href="mailto:minniti@ipmc.cnrs.fr">minniti@ipmc.cnrs.fr</a> or via 
        <a href="https://github.com/Jumitti/LabMaster/issues">Issues</a>.
    </div>
    ''',
    unsafe_allow_html=True
)
st.markdown("")
