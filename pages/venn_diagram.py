import csv
import io
from io import BytesIO
from itertools import combinations
from zipfile import ZipFile

import chardet
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import streamlit as st
from venn import venn, pseudovenn
import altair_upset as au
import altair as alt

from utils.page_config import page_config
import os


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
def download_png(graph_type="venn"):
    buffer_png = BytesIO()
    if graph_type == "venn":
        plt.savefig(buffer_png, format="png", bbox_inches='tight')
    elif graph_type == "upset":
        chart.save(buffer_png, format="png")
    buffer_png.seek(0)

    return buffer_png


# For download SVG Venn
def download_svg(graph_type="venn"):
    buffer_svg = BytesIO()
    if graph_type == "venn":
        plt.savefig(buffer_svg, format="svg", bbox_inches='tight')
    elif graph_type == "upset":
        chart.save(buffer_svg, format="svg")
    buffer_svg.seek(0)

    return buffer_svg


# Example files
csv_file = os.path.join(os.path.dirname(__file__), '../example', 'example_venn_diagram.csv')
xlsx_file = os.path.join(os.path.dirname(__file__), '../example', 'example_venn_diagram.xlsx')

# Setting of Venn configurations
fmt_options = {"Number": "{size}",
               "Percentage": "{percentage:.1f}%",
               "Logic": "{logic}"}

cmap_options = {'Accent': 'Accent',
                'BRG': 'brg',
                "Civids": 'cividis',
                'CMRmap': 'CMRmap',
                "CoolWarm": "coolwarm",
                'CubeHelix': 'cubehelix',
                'Dark2': 'Dark2',
                "Default": 'hsv',
                'Flag': 'flag',
                'Gist Earth': 'gist_earth',
                'Gist Ncar': 'gist_ncar',
                'Gist Rainbow': 'gist_rainbow',
                'Gist Stern': 'gist_stern',
                'GnuPlot': 'gnuplot',
                'GnuPlot2': 'gnuplot2',
                "Inferno": 'inferno',
                'Jet': 'jet',
                "Magma": 'magma',
                'Nipy Spectral': 'nipy_spectral',
                'Ocean': 'ocean',
                'Pastel2': 'Pastel2',
                'Paired': 'Paired',
                "Plasma": 'plasma',
                'Prism': 'prism',
                'Rainbow': 'rainbow',
                'Set1': 'Set1',
                'Set2': 'Set2',
                'Set3': 'Set3',
                "Spectral": "Spectral",
                'Tab10': 'tab10',
                'Tab20': 'tab20',
                'Tab20b': 'tab20b',
                'Tab20c': 'tab20c',
                'Terrain': 'terrain',
                'Turbo': 'turbo',
                "Twilight": 'twilight',
                "Twilight Shifted": 'twilight_shifted',
                "Viridis": 'viridis',
                }

legend_loc_options = {'Best': 'best',
                      'Upper Right': 'upper right',
                      'Upper Left': 'upper left',
                      'Upper Center': 'upper center',
                      'Lower Right': 'lower right',
                      'Lower Left': 'lower left',
                      'Lower Center': 'lower center',
                      'Right': 'right',
                      'Center Right': 'center right',
                      'Center Left': 'center left',
                      'Center': 'center'}

# Setting of UpSet plot configurations
sorted_by_option = {"Frequency": "frequency",
                    "Degree": "degree"}

sorted_order_option = {"Descending": "descending",
                       "Ascending": "ascending"}

# Page config
page_config()

# Main page
st.title('⭕ VennLit V2')

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

            with col2:
                # Lists selection
                st.subheader('📌 Lists selection')
                items_occurrence = {per_list: set(df[per_list].dropna()) for per_list in lists}
                selection_lists = st.multiselect('Lists selection', lists, default=lists[:2],
                                                 placeholder="Choose lists", disabled=False,
                                                 label_visibility='collapsed')
                num_sets = len(selection_lists)
                selected_lists = selection_lists[:num_sets]
    except Exception as e:
        with col2:
            st.warning(f"It appears that there is an error with one or more values in your lists..."
                       "Please check your data. Otherwise, convert your file to .csv with the ';' deliminator.\n\n"
                       "If this does not resolve the problems, contact me by email (minnitijulien06@gmail.com ; minniti@ipmc.cnrs.fr) or submit a [GitHub Issue](https://github.com/Jumitti/vennlit_v2/issues).\n\n"
                       "Error information:"
                       f"{e}", icon='🚨')

    with col1:
        # Credits section
        with st.expander("✒️Credits", expanded=False):
            st.write("Original app by [@professordata](https://github.com/dataprofessor/vennlit)")
            st.write(
                "Venn diagram with [@tctianchi](https://github.com/tctianchi/pyvenn) and [@LankyCyril](https://github.com/LankyCyril/pyvenn)")
            st.write(
                "Inspired by [InteractiVenn](http://www.interactivenn.net/) (DOI:[10.1186/s12859-015-0611-3](http://doi.org/10.1186/s12859-015-0611-3))")
            st.write("VennLit V2 rebuild and up-to-date by [@Jumitti](https://github.com/Jumitti/vennlit_v2)")

try:
    plt.figure(figsize=(8, 8))
    if len(selection_lists) > 1:
        with col3:
            st.subheader("")
            venn_data = download_venn_data(selected_lists)
            st.download_button(label="💾 Download Venn data",
                               data=venn_data,
                               file_name=f'venn_data{"".join("_" + selected_list for selected_list in selection_lists)}.zip',
                               mime="application/zip", )

    if 1 < len(selection_lists) <= 6:  # Venn diagram 2 to 6 comparisons
        with col3:
            st.divider()
            # Settings
            st.subheader('⚙️Venn diagram settings')
            fmt = st.radio(
                "**Number format:**",
                list(fmt_options.keys()),
                index=0,
                horizontal=True, key='venn_fmt')
            venn_format = fmt_options[fmt]

            cmap = st.selectbox(
                "**Colors:**",
                list(cmap_options.keys()),
                index=7, key='venn_cmap')
            cmap_format = cmap_options[cmap]

            font_size = st.slider("**Font size:**", min_value=5, max_value=20, value=10, step=1, key='venn_font_size',
                                  help=None)

            fig_size = st.slider("**Venn size**:", min_value=5, max_value=20, value=10, step=1, key='venn_fig_size',
                                 help=None)

            legend_loc = st.selectbox(
                "**Legend position:**",
                list(legend_loc_options.keys()),
                index=0, key='venn_legend_loc')
            legend_loc_format = legend_loc_options[legend_loc]

        with col2:
            # Venn diagram
            st.subheader('Venn diagram')
            dataset_dict = {name: set(items_occurrence[name]) for name in selected_lists}
            venn(dataset_dict, fmt=venn_format, cmap=cmap_format, fontsize=font_size, legend_loc=legend_loc_format,
                 figsize=(fig_size, fig_size))
            st.pyplot(plt)

        with col3:
            # Download PNG and SVG
            buffer_png = download_png()
            st.download_button(
                label="💾 Download Venn diagram (.png)",
                data=buffer_png,
                file_name=f'venn{"".join("_" + selected_list for selected_list in selection_lists)}.png',
                mime='image/png',
            )

            buffer_svg = download_svg()
            st.download_button(
                label="💾 Download Venn diagram (.svg)",
                data=buffer_svg,
                file_name=f'venn{"".join("_" + selected_list for selected_list in selection_lists)}.svg',
                mime='image/svg+xml',
            )
            st.write(
                'Try opening the .svg diagram using [Inkscape](https://inkscape.org/) to move shapes, resize, change font, colors and more.')

    if len(selection_lists) == 6:  # Pseudo-Venn for 6 comparison
        with col3:
            st.divider()
            # Pseudo-Venn settings
            st.subheader('⚙️Pseudo-Venn diagram settings',
                         help='Six-set true Venn diagrams are somewhat unwieldy, and not all intersections are usually of interest.\n\n'
                              'If you wish to display information about elements in hidden intersections,'
                              'uncheck the option **hidden intersections** below.\n\n'
                              'Some intersections are not present, but the most commonly wanted are.')
            fmt = st.radio(
                "**Number format:**",
                list(fmt_options.keys()),
                index=0,
                horizontal=True, key='pseudovenn_fmt')
            venn_format = fmt_options[fmt]

            cmap = st.selectbox(
                "**Colors:**",
                list(cmap_options.keys()),
                index=7, key='pseudovenn_cmap')
            cmap_format = cmap_options[cmap]

            font_size = st.slider("**Font size:**", min_value=5, max_value=20, value=10, step=1,
                                  key='pseudovenn_font_size',
                                  help=None)

            fig_size = st.slider("**Pseudo-Venn size:**", min_value=5, max_value=20, value=10, step=1,
                                 key='pseudovenn_fig_size',
                                 help=None)

            legend_loc = st.selectbox(
                "**Legend position:**",
                list(legend_loc_options.keys()),
                index=0, key='pseudovenn_legend_loc')
            legend_loc_format = legend_loc_options[legend_loc]

            hint_hidden_format = st.checkbox('**Hidden intersections**', value=1,
                                             help='Six-set true Venn diagrams are somewhat unwieldy, and not all intersections are usually of interest.\n\n'
                                                  'If you wish to display information about elements in hidden intersections,'
                                                  'uncheck the option **hidden intersections**.\n\n'
                                                  'Some intersections are not present, but the most commonly wanted are.')

        with col2:
            # Pseudo-Venn diagram
            st.subheader('Pseudo-Venn diagram')
            dataset_dict = {name: set(items_occurrence[name]) for name in selected_lists}
            pseudovenn(dataset_dict, fmt=venn_format, cmap=cmap_format, fontsize=font_size,
                       legend_loc=legend_loc_format,
                       figsize=(fig_size, fig_size),
                       hint_hidden=False if hint_hidden_format else True)
            st.pyplot(plt)

        with col3:
            # Download PNG and SVG
            buffer_png = download_png()
            st.download_button(
                label="💾 Download Pseudo-Venn diagram (.png)",
                data=buffer_png,
                file_name=f'pseudovenn{"".join("_" + selected_list for selected_list in selection_lists)}.png',
                mime='image/png',
            )

            buffer_svg = download_svg()
            st.download_button(
                label="💾 Download Pseudo-Venn diagram (.svg)",
                data=buffer_svg,
                file_name=f'pseudovenn{"".join("_" + selected_list for selected_list in selection_lists)}.svg',
                mime='image/svg+xml',
            )
            st.write(
                'Try opening the .svg diagram using [Inkscape](https://inkscape.org/) to move shapes, resize, change font, colors and more.')

    # Upset plot
    with col3:
        st.divider()
        # Settings
        st.subheader('⚙️UpSet plot settings')
        sorted_by = st.radio(
            "**Number format:**",
            list(sorted_by_option.keys()),
            index=0,
            horizontal=True)
        sorted_by = sorted_by_option[sorted_by]

        sorted_order = st.radio(
            "**Number format:**",
            list(sorted_order_option.keys()),
            index=0,
            horizontal=True)
        sorted_order = sorted_order_option[sorted_order]

        color = st.color_picker("**Main/highlight color:**", "#215FD2")

        cmap = st.selectbox(
            "**Sets colors:**",
            list(cmap_options.keys()),
            index=7, key='upset_cmap')
        cmap_format = cmap_options[cmap]
        cmap = plt.get_cmap(cmap_format)
        colors_hex = [mcolors.to_hex(cmap(i / (len(selection_lists) - 1))) for i in range(len(selection_lists))]

    with col2:
        st.subheader('UpSet plot')
        data = {name: [] for name in selected_lists}
        all_items = set.union(*[set(items_occurrence[k]) for k in selected_lists])

        for item in all_items:
            for key in selected_lists:
                data[key].append(1 if item in items_occurrence[key] else 0)

        data_upset = pd.DataFrame(data)
        chart = au.UpSetAltair(data=data_upset, sets=selected_lists, sort_by=sorted_by, sort_order=sorted_order,
                               highlight_color=color, color_range=colors_hex)
        buffer_png = download_png(graph_type="upset")

        st.altair_chart(chart, use_container_width=False)

    col3.download_button(
        label="💾 Download UpSet plot (.png)",
        data=buffer_png,
        file_name=f'upsetplot{"".join("_" + selected_list for selected_list in selection_lists)}.png',
        mime="image/png"
    )

except Exception as e:
    with col2:
        st.warning(f"It appears that there is an error with one or more values in your lists..."
                   "Please check your data. Otherwise, convert your file to .csv with the ';' deliminator.\n\n"
                   "If this does not resolve the problems, contact me by email (minnitijulien06@gmail.com ; minniti@ipmc.cnrs.fr) or submit a [GitHub Issue](https://github.com/Jumitti/LabMaster/issues).\n\n"
                   "Error information:"
                   f"{e}", icon='🚨')
