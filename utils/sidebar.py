import streamlit as st
import platform
import math
import os

import pandas as pd

local_test = platform.processor()


def sidebar():
    # Main page
    img_path = os.path.join(os.path.dirname(__file__), './img', 'labmaster_logo.png')

    st.logo(img/labmaster_logo.png)
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
    st.sidebar.page_link("pages/venn_diagram.py", label="Venn diagram")

    # MIT licence
    st.sidebar.divider()
    st.sidebar.link_button('Under MIT licence', 'https://github.com/Jumitti/BlotMaster/blob/main/LICENSE')