import csv
import io
import os
import warnings
import time
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


# Example files
csv_file = os.path.join(os.path.dirname(__file__), '../example', 'example_cohen_s_d.csv')
xlsx_file = os.path.join(os.path.dirname(__file__), '../example', 'example_cohen_s_d.xlsx')

# Page config
page_config()

# Main page
st.title("Effect Size Calculator")

df = []
selection_lists = []

col1, col2, col3 = st.columns([0.8, 1.4, 0.8])

with col1:
    # Example section
    st.subheader("📎 Example and Hints")

    st.link_button("Help ?!", 'https://pingouin-stats.org/build/html/generated/pingouin.compute_effsize.html#pingouin.compute_effsize')

    demo = st.checkbox("**Try example**", value=1)
    if demo:  # Demo mode
        with col2:
            st.write('You are by default in **demo** mode.\n'
                     'You can play with **Demo** or disable **Try example** on the left **📎 Example** section.\n'
                     'You can also click on **[Help](https://pingouin-stats.org/build/html/generated/pingouin.compute_effsize.html#pingouin.compute_effsize)**.')
        snif_delimiter = detect_delimiter(csv_file)
        df = pd.read_csv(csv_file, delimiter=snif_delimiter)

    with st.expander("**.csv and .xlsx templates**", expanded=False):
        with open(csv_file, "rb") as file:  # Download .csv template
            st.download_button(
                label="Download example.csv",
                data=file,
                file_name="example_cohen_s_d.csv",
                mime="text/csv")
        with open(xlsx_file, "rb") as file:  # Download .xlsx template
            st.download_button(
                label="Download example.xlsx",
                data=file,
                file_name="example_cohen_s_d.xlsx",
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
            st.write('You can play with **Try example** on the left **📎 Example** section.\n'
                     'You can also click on **[Help](https://pingouin-stats.org/build/html/generated/pingouin.compute_effsize.html#pingouin.compute_effsize)**.')
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

            with col2:
                col2a, col2b = st.columns([2, 1])
                if len(selected_lists) > 1:
                    col3.subheader("**Effect size settings**", help=f"Please, see [HELP](https://pingouin-stats.org/build/html/generated/pingouin.compute_effsize.html#pingouin.compute_effsize)")

                    n = len(selected_lists)
                    results_matrix = np.zeros((n, n))

                    col3a, col3b = col3.columns(2)
                    col3a.markdown("Effect size type")
                    methode = col3a.radio("Effect size type",
                                          ("Cohen's d (unpaired)", "Cohen's d' (paired)", "Hedges g (unpaired)", "Hedges g (paired)"),
                                          label_visibility="collapsed")

                    col3b.markdown("Reverse comparison")
                    invert = col3b.toggle("Reverse comparison", value=False)
                    if not invert:
                        col3b.write(f"**Example**: {selected_lists[0]} vs. {selected_lists[1]}")
                    else:
                        col3b.write(f"**Example**: {selected_lists[1]} vs. {selected_lists[0]}")

                    with warnings.catch_warnings(record=True) as w:
                        warnings.simplefilter("always")

                        different_group_size = []
                        short_n = 0

                        for i in range(n):
                            for j in range(i + 1, n):
                                group1 = df[selected_lists[i]].dropna()
                                group2 = df[selected_lists[j]].dropna()

                                if len(group1) < 20 or len(group2) < 20:
                                    short_n += 1

                                if invert:
                                    a, b = group1, group2
                                else:
                                    a, b = group2, group1

                                if methode == "Cohen's d (unpaired)":
                                    d_result = pg.compute_effsize(a, b, eftype='cohen')
                                elif methode == "Cohen's d' (paired)":
                                    d_result = pg.compute_effsize(a, b, paired=True, eftype='cohen')
                                elif methode == "Hedges g (unpaired)":
                                    d_result = pg.compute_effsize(a, b, eftype='hedges')
                                elif methode == "Hedges g (paired)":
                                    d_result = pg.compute_effsize(a, b, paired=True, eftype='hedges')

                                if len(w) > 0:
                                    for warning in w:
                                        if "Switching to paired == False" in str(warning.message):
                                            different_group_size.append((selected_lists[i], selected_lists[j]))
                                            if methode == "Cohen's d' (paired)":
                                                d_result = pg.compute_effsize(a, b, eftype='cohen')
                                            elif methode == "Hedges g (paired)":
                                                d_result = pg.compute_effsize(a, b, eftype='hedges')

                                results_matrix[i, j] = d_result
                                results_matrix[j, i] = d_result

                    if short_n > 0 and methode in ["Cohen's d (unpaired)", "Cohen's d' (paired)"]:
                        col3.warning("The Cohen’s d is a biased estimate of the population effect size, especially for "
                                     "small samples (n < 20). It is often preferable to use the corrected Hedges g instead")

                    if len(different_group_size) > 0:
                        col3.warning(f"Please note some groups are not the same size. Analysis was automatically performed as unpaired groups.")

                    col2a.subheader("Effect size")
                    results_df = pd.DataFrame(results_matrix, columns=selected_lists, index=selected_lists)
                    results_df = results_df.sort_index(axis=0).sort_index(axis=1)
                    col2a.dataframe(results_df.style.background_gradient(cmap="coolwarm", axis=None))

                    mean = [df[col].mean() for col in selected_lists]
                    sd = [df[col].std() for col in selected_lists]
                    n = [len(df[col].dropna()) for col in selected_lists]

                    descriptive_df = pd.DataFrame({
                        'Groups': selected_lists,
                        'Mean': mean,
                        'SD': sd,
                        'n': n
                    })

                    col2b.subheader("Desc. Statistics")
                    col2b.dataframe(descriptive_df.sort_values(by="Groups"), hide_index=True)

                    data = pd.DataFrame({
                        'Groups': selected_lists,
                        'Mean': mean,
                        'SD': sd
                    })

                    barplot = alt.Chart(data).mark_bar().encode(
                        x='Groups',
                        y='Mean',
                        color='Groups'
                    ).properties().interactive()

                    error_bar = barplot.mark_errorbar().encode(
                        y=alt.Y('Mean:Q'),
                        yError='SD'
                    )

                    st.altair_chart(barplot + error_bar, use_container_width=True)

                else:
                    st.info("Please select at least two groups to compare.")

            with col3:
                st.divider()
                if len(selected_lists) > 2:
                    st.subheader("Eta-square settings", help=f"Please, see [HELP](https://pingouin-stats.org/build/html/generated/pingouin.anova.html)")
                    # Champ de saisie pour le nombre de listes à comparer
                    num_groups = st.number_input("Number of groups to compare (minimum 3)", min_value=3, step=1)
                    detailed = st.toggle("Detailed", value=False)

                    # Liste de couleurs visibles et pâles
                    colors = [
                        'rgba(255, 99, 132, 0.3)',  # Rose pâle
                        'rgba(54, 162, 235, 0.3)',  # Bleu pâle
                        'rgba(255, 206, 86, 0.3)',  # Jaune pâle
                        'rgba(75, 192, 192, 0.3)',  # Vert pâle
                        'rgba(153, 102, 255, 0.3)'  # Violet pâle
                    ]

                    if len(selected_lists) >= num_groups:
                        # Prendre les combinaisons de groupes selon le nombre spécifié
                        group_combinations = list(combinations(selected_lists, num_groups))

                        # Initialiser un DataFrame pour stocker les résultats
                        results_list = []

                        for index, group_set in enumerate(group_combinations):
                            # Effectuer l'ANOVA pour les groupes choisis
                            long_df = df[list(group_set)].melt(var_name="Groups", value_name="Values").dropna()
                            anova_result = pg.anova(data=long_df, dv="Values", between="Groups", detailed=detailed, effsize="n2")

                            # Extraire l'eta-squared directement des résultats de l'ANOVA
                            eta_squared_value = anova_result['n2'][
                                0]  # La première ligne contient l'eta-squared pour l'ANOVA

                            # Ajouter une colonne pour eta-squared dans le résultat ANOVA
                            anova_result['Eta-Squared'] = [eta_squared_value] + [''] * (
                                        len(anova_result) - 1)  # Ajouter à la première ligne, le reste vide

                            # Enlever les index non désirés
                            anova_result.reset_index(drop=True, inplace=True)

                            # Ajouter les résultats dans une liste
                            results_list.append({
                                "Groups": ', '.join(group_set),
                                "ANOVA Result": anova_result,
                                "Color": colors[index % len(colors)]  # Associer une couleur à chaque groupe
                            })

                        # Créer un DataFrame pour tous les résultats
                        all_results = pd.concat([result['ANOVA Result'] for result in results_list],
                                                keys=[result['Groups'] for result in results_list], names=['Groups'])

                        # Réinitialiser l'index pour enlever les index indésirables
                        all_results.reset_index(level=0,
                                                inplace=True)  # Cela déplace le niveau 'Groups' au niveau de la colonne

                        # Supprimer la colonne d'index qui contient 0 et 1
                        all_results.index = range(len(all_results))  # Réinitialiser l'index

                        # Afficher les résultats consolidés avec couleurs
                        col2.subheader("ANOVA and Eta-square results", help=f"Please, see [HELP](https://pingouin-stats.org/build/html/generated/pingouin.anova.html)")

                        # Créer une fonction de coloration
                        def highlight_groups(row):
                            # Utiliser l'index pour assigner la couleur
                            color = colors[
                                results_list.index(next(r for r in results_list if r['Groups'] == row['Groups'])) % len(
                                    colors)]
                            return [f'background-color: {color}' for _ in row]


                        # Appliquer le style sans afficher la colonne "Color"
                        styled_results = all_results.style.apply(highlight_groups, axis=1)

                        # Afficher le DataFrame
                        col2.dataframe(styled_results, hide_index=True, use_container_width=True)

                    else:
                        col3.info("Please select at least as many groups as the number you wish to compare.")
    except Exception as e:
        with col2:
            st.warning(f"It appears that there is an error with one or more values in your lists..."
                       "Please check your data. Otherwise, convert your file to .csv with the ';' deliminator.\n\n"
                       "If this does not resolve the problems, contact me by email (minnitijulien06@gmail.com ; minniti@ipmc.cnrs.fr) or submit a [GitHub Issue](https://github.com/Jumitti/LabMaster/issues).\n\n"
                       "Error information:"
                       f"{e}", icon='🚨')
