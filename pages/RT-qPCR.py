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

# MIT licence
st.sidebar.divider()
st.sidebar.link_button('Under MIT licence', 'https://github.com/Jumitti/BlotMaster/blob/main/LICENSE')