import datetime
import io
import json

from Bio import SeqIO
from io import StringIO

import altair as alt
import pandas as pd
import requests
import streamlit as st
from stqdm import stqdm
from pages.design_primer_API import NCBIdna, Primer3

from utils.page_config import page_config


ucsc_species = {
    "Homo sapiens": {'org': 'Human', 'db': 'hg38', 'wp_target': ['genome', 'hg38KgSeqV48']},
    "Mus musculus": {'org': 'Mouse', 'db': 'mm39', 'wp_target': ['genome', 'mm39KgSeqVM37']}
}


def is_valid_dna(seq):
    return all(base in "ATGCatgc" for base in seq)


def reverse_complement(seq):
    complement = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(complement)[::-1]


# Page config
page_config()

if st.button("You can also design your own primers by clicking here ☺"):
    st.switch_page("pages/design_primer.py")

if "primers" not in st.session_state:
    st.session_state.primers = []

st.header("🧬 Primer Submission")

with st.form("primer_form"):
    col1, col2, col3 = st.columns(3)
    with col1:
        forward = st.text_input("Forward Primer (5'→3')")
    with col2:
        reverse = st.text_input("Reverse Primer (5'→3')")
        flip_reverse = st.toggle("↩️ Flip reverse primer (reverse complement)")
    with col3:
        species = st.selectbox("Species", ["Homo sapiens", "Mus musculus", "Other/Unknown"])

    submitted = st.form_submit_button("➕ Add Primer")

    if submitted:
        if not forward or not reverse:
            st.error("❌ Both primers must be filled.")
        elif not is_valid_dna(forward) or not is_valid_dna(reverse):
            st.error("❌ Only A, T, G, C characters are allowed.")
        else:
            if flip_reverse:
                reverse = reverse_complement(reverse)
            st.session_state.primers.append({
                "Forward": forward.upper(),
                "Reverse": reverse.upper(),
                "Species": species
            })
            st.success("✅ Primer added!")

if st.session_state.primers:
    df = pd.DataFrame(st.session_state.primers)

    df["❌ Delete"] = False

    column_config = {
        "Species": st.column_config.SelectboxColumn(
            "Species", options=["Homo sapiens", "Mus musculus", "Other/Unknown"], required=True
        ),
        "Forward": st.column_config.TextColumn("Forward", disabled=True),
        "Reverse": st.column_config.TextColumn("Reverse", disabled=True),
        "❌ Delete": st.column_config.CheckboxColumn("Delete", help="Check to remove this primer")
    }

    st.write("### 🧬 Submitted Primers")
    edited_df = st.data_editor(
        df,
        column_config=column_config,
        use_container_width=True,
        hide_index=False,
        key="primer_editor"
    )

    col_a, col_b = st.columns([3, 1])
    with col_a:
        if st.button("🎯🗑️ Clear Selected Primers"):
            new_data = edited_df[edited_df["❌ Delete"] == False].drop(columns="❌ Delete")
            st.session_state.primers = new_data.to_dict(orient="records")
            st.success("Selected primers removed.")
            st.rerun()
    with col_b:
        if st.button("🗑️ Clear All Primers"):
            st.session_state.primers.clear()
            st.success("All primers cleared.")
            st.rerun()

    if st.button("Run"):
        if species in ucsc_species.keys() and ucsc_validation is True:
            org = ucsc_species[species]["org"]
            db = ucsc_species[species]["db"]
            wp_targets = ucsc_species[species]["wp_target"]
            validation_relative, validation_absolute, sequence_relative, sequence_absolute = Primer3.fetch_ucsc_pcr_results(
                species, org, db, wp_targets, left_seq, right_seq,
                amplicon_size_abs, PRIMER_PRODUCT_SIZE_RANGE[1])

else:
    st.info("No primers submitted yet.")


if st.button("You can also design your own primers by clicking here ☺", key="design_primer"):
    st.switch_page("pages/design_primer.py")