import csv
import io
import os
from io import BytesIO
from itertools import combinations
from zipfile import ZipFile

import altair as alt
import chardet
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pingouin as pg
import streamlit as st

from utils.page_config import page_config


@st.cache_data(ttl=3600)
def load_sheet_names(file):
    df = pd.ExcelFile(file)
    sheet_names = df.sheet_names
    return sheet_names


@st.cache_data(ttl=3600)
def load_data(file, selected_sheet=None):
    if file.type == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        df = pd.read_excel(file, sheet_name=selected_sheet)
    else:
        file_content = file.read()
        delimiter, encoding = detect_delimiter(file_content)
        df = pd.read_csv(io.StringIO(file_content.decode(encoding)), delimiter=delimiter)
    return df


@st.cache_data(ttl=3600)
def detect_delimiter(file_path):
    try:
        with open(file_path, 'rb') as f:
            result = chardet.detect(f.read())
        encoding = result['encoding']
        with open(file_path, 'r', encoding=encoding) as f:
            sample = f.read(4096)
            dialect = csv.Sniffer().sniff(sample)
        return dialect.delimiter
    except Exception as e:
        result = chardet.detect(file_path)
        encoding = result['encoding']
        sample = file_path[:4096]
        sample = sample.decode(encoding)
        dialect = csv.Sniffer().sniff(sample)

        return dialect.delimiter, encoding


# Analyse of selected list to generate multidimensional Venn files
@st.cache_data(ttl=3600)
def download_venn_data(lists):
    items_occurrence = {per_list: set(df[per_list].dropna()) for per_list in lists}
    zip_buffer = BytesIO()

    with ZipFile(zip_buffer, 'a') as zip_file:
        for current_list, items_current_list in items_occurrence.items():
            exclusive_items = items_current_list.copy()

            for other_list, items_other_list in items_occurrence.items():
                if other_list != current_list:
                    exclusive_items -= items_other_list

            file_content = "\n".join(map(str, exclusive_items))
            zip_file.writestr(f"1_{current_list}.txt", file_content)

        for combination_size in range(2, len(lists) + 1):
            for lists_combination in combinations(lists, combination_size):
                shared_items = set.intersection(*(items_occurrence[item] for item in lists_combination))
                items_exclusive_to_combination = shared_items.copy()

                for other_list, items_other_list in items_occurrence.items():
                    if other_list not in lists_combination:
                        items_exclusive_to_combination -= items_other_list

                file_content = "\n".join(items_exclusive_to_combination)
                file_name = f"{len(lists_combination)}_{'_'.join(sorted(lists_combination))}.txt"
                zip_file.writestr(file_name, file_content)

    venn_data = zip_buffer.getvalue()
    return venn_data


# For download PNG Venn
def download_png():
    buffer_png = BytesIO()
    plt.savefig(buffer_png, format="png", bbox_inches='tight')
    buffer_png.seek(0)

    return buffer_png


# For download SVG Venn
def download_svg():
    buffer_svg = BytesIO()
    plt.savefig(buffer_svg, format="svg", bbox_inches='tight')
    buffer_svg.seek(0)

    return buffer_svg


# Example files
csv_file = os.path.join(os.path.dirname(__file__), '../example', 'example_cohen_s_d.csv')
xlsx_file = os.path.join(os.path.dirname(__file__), '../example', 'example_cohen_s_d.xlsx')

# Page config
page_config()

# Titre de l'application
st.title("Comparaison des groupes avec Cohen's d")

df = []
selection_lists = []

col1, col2, col3 = st.columns([0.8, 1.4, 0.8])

with col1:
    # Example section
    st.subheader("📎 Example and Hints")

    st.link_button("Help", 'https://jumitti.notion.site/jumitti/VennLit-V2-e20a373a9c6f4c1390e72a7953ffcb0c')

    demo = st.checkbox("**Try example**", value=1)
    if demo:  # Demo mode
        with col2:
            st.subheader('Welcome to VennLit V2 😊')
            st.write('You are by default in **demo** mode.\n'
                     'You can play with VennLit V2 or disable **Try example** on the left **📎 Example** section.\n'
                     'You can also click on **[Help](https://jumitti.notion.site/jumitti/VennLit-V2-e20a373a9c6f4c1390e72a7953ffcb0c)**.')
        snif_delimiter = detect_delimiter(csv_file)
        df = pd.read_csv(csv_file, delimiter=snif_delimiter)

    with st.expander("**.csv and .xlsx templates**", expanded=False):
        with open(csv_file, "rb") as file:  # Download .csv template
            st.download_button(
                label="Download example.csv",
                data=file,
                file_name="example_venn_diagram.csv",
                mime="text/csv")
        with open(xlsx_file, "rb") as file:  # Download .xlsx template
            st.download_button(
                label="Download example.xlsx",
                data=file,
                file_name="example_venn_diagram.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    # Upload data section
    st.subheader("💽 Upload data")

    uploaded_files = st.file_uploader("**Upload one or more .xlsx .csv files**", type=["csv", "xlsx"],
                                      accept_multiple_files=True)
    if len(uploaded_files) > 0:
        dfs = []
        for file in uploaded_files:  # Is .csv or .xlsx
            if file.type == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
                sheet_names = load_sheet_names(file)
                selected_sheet = st.selectbox(f"Select sheet for **{file.name}**:", sheet_names)
                df = load_data(file, selected_sheet=selected_sheet)
            else:
                df = load_data(file)

            dfs.append(df)

        all_columns = [col for df in dfs for col in df.columns]
        duplicate_columns = [col for col in set(all_columns) if all_columns.count(col) > 1]
        if duplicate_columns:  # Some lists with same name ?
            st.warning(f"Some lists have the same name: {', '.join(duplicate_columns)}")
            filtered_dfs = []
            included_columns = set()
            for df in dfs:
                filtered_columns = [col for col in df.columns if col not in duplicate_columns]
                included_columns.update(filtered_columns)
                filtered_df = df[filtered_columns]
                filtered_dfs.append(filtered_df)
            df = pd.concat(filtered_dfs, axis=1)
        else:
            df = pd.concat(dfs, axis=1)
    elif len(uploaded_files) == 0 and not demo:
        st.cache_data.clear()
        with col2:
            st.subheader('Welcome to VennLit V2 😊')
            st.write('You can play with VennLit V2 or enable **Try example** on the left **📎 Example** section.\n'
                     'You can also click on **[Help](https://jumitti.notion.site/jumitti/VennLit-V2-e20a373a9c6f4c1390e72a7953ffcb0c)**.')
    else:
        st.cache_data.clear()

    try:
        if len(df) > 0:
            # Lists section
            st.subheader("🧮 Lists")
            st.dataframe(df, hide_index=True)
            lists = df.columns.tolist()

            with col3:
                # Lists selection
                st.subheader('📌 Lists selection')
                items_occurrence = {per_list: set(df[per_list].dropna()) for per_list in lists}
                selection_lists = st.multiselect('Lists selection', lists, default=lists[:2],
                                                 placeholder="Choose lists", disabled=False,
                                                 label_visibility='collapsed')
                num_sets = len(selection_lists)
                selected_lists = selection_lists[:num_sets]

                # Radio button pour sélectionner le type d'effet de Cohen
                methode = st.radio("Effect size type",
                                   ("Cohen's d (unpaired)", "Cohen's d (paired)"))

            with col2:
                if len(selected_lists) >= 2:
                    # Initialiser un tableau pour stocker les résultats
                    n = len(selected_lists)
                    results_matrix = np.zeros((n, n))

                    # Calculer Cohen's d ou d' pour chaque paire de colonnes
                    for i in range(n):
                        for j in range(i + 1, n):
                            groupe1 = df[selected_lists[i]].dropna()
                            groupe2 = df[selected_lists[j]].dropna()

                            if methode == "Cohen's d (unpaired)":
                                d_result = pg.compute_effsize(groupe1, groupe2, paired=False, eftype='cohen')
                                print("Cohen's d (unpaired)")
                            else:
                                d_result = pg.compute_effsize(groupe1, groupe2, paired=True, eftype='cohen')

                            # Stocker les résultats dans la matrice
                            results_matrix[i, j] = d_result
                            results_matrix[j, i] = d_result

                    # Afficher la matrice des tailles d'effet
                    st.subheader("Échiquier des tailles d'effet (Cohen's d)")
                    results_df = pd.DataFrame(results_matrix, columns=selected_lists, index=selected_lists)
                    st.dataframe(results_df.style.background_gradient(cmap="coolwarm", axis=None))

                    # Calculer les moyennes et erreurs standards pour chaque groupe
                    moyennes = [df[col].mean() for col in selected_lists]
                    se = [df[col].std() / (len(df[col].dropna()) ** 0.5) for col in selected_lists]

                    # Créer un DataFrame pour Altair (visualisation des moyennes)
                    data = pd.DataFrame({
                        'Groupe': selected_lists,
                        'Moyenne': moyennes,
                        'SE': se
                    })

                    # Créer un barplot avec Altair
                    barplot = alt.Chart(data).mark_bar().encode(
                        x='Groupe',
                        y='Moyenne',
                        color='Groupe'
                    ).properties(
                        title=f"Comparaison des groupes avec Cohen's d ou d'"
                    )

                    # Ajouter des barres d'erreur
                    erreur = barplot.mark_errorbar().encode(
                        y=alt.Y('Moyenne:Q'),
                        yError='SE'
                    )

                    # Afficher le graphique
                    st.altair_chart(barplot + erreur, use_container_width=True)

                else:
                    st.info("Veuillez sélectionner au moins deux colonnes à comparer.")
    except Exception as e:
        with col2:
            st.warning(f"It appears that there is an error with one or more values in your lists..."
                       "Please check your data. Otherwise, convert your file to .csv with the ';' deliminator.\n\n"
                       "If this does not resolve the problems, contact me by email (minnitijulien06@gmail.com ; minniti@ipmc.cnrs.fr) or submit a [GitHub Issue](https://github.com/Jumitti/vennlit_v2/issues).\n\n"
                       "Error information:"
                       f"{e}", icon='🚨')
