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
from pages.design_primer_API import NCBIdna

from utils.page_config import page_config


def graphique(exons, primers, normalization=False):
    if normalization is True:
        adjusted_exons = []
        intron_cumulative_length = 0
        for i, (start, end) in enumerate(exons):
            if i > 0:
                intron_cumulative_length += exons[i][0] - exons[i - 1][1]

            adjusted_start = start - intron_cumulative_length
            adjusted_end = end - intron_cumulative_length
            adjusted_exons.append((adjusted_start, adjusted_end))

    exons_data = pd.DataFrame(adjusted_exons if normalization is True else exons, columns=['start', 'end'])
    exons_data['exon_number'] = exons_data.index + 1

    primers_data = [{'pair': f'Pair {i + 1}', 'type': 'Left Primer',
                     'sequence': p['left_primer']['sequence'],
                     'length': p['left_primer']['length'],
                     'position': f"{p['left_primer']['position'][0]}-{p['left_primer']['position'][1]}",
                     'position_start': p['left_primer']['position'][0],
                     'position_end': p['left_primer']['position'][1],
                     'position_abs': f"{p['left_primer']['position_abs'][0]}-{p['left_primer']['position_abs'][1]}",
                     'position_start_abs': p['left_primer']['position_abs'][0],
                     'position_end_abs': p['left_primer']['position_abs'][1],
                     'tm': p['left_primer']['tm'],
                     'self_complementarity': p['left_primer']['self_complementarity'],
                     "self_3prime_complementarity": p['left_primer']['self_3prime_complementarity'],
                     'amplicon_size': p['amplicon_size'],
                     'amplicon_size_abs': p['amplicon_size_abs'],
                     'label': f'Left {i + 1}', 'y': (i + 1) * 10}
                    for i, p in enumerate(primers)] + [
        {'pair': f'Pair {i + 1}', 'type': 'Right Primer',
         'sequence': p['right_primer']['sequence'],
         'length': p['right_primer']['length'],
         'position': f"{p['right_primer']['position'][0]}-{p['right_primer']['position'][1]}",
         'position_start': p['right_primer']['position'][0],
         'position_end': p['right_primer']['position'][1],
         'position_abs': f"{p['right_primer']['position_abs'][0]}-{p['right_primer']['position_abs'][1]}",
         'position_start_abs': p['right_primer']['position_abs'][0],
         'position_end_abs': p['right_primer']['position_abs'][1],
         'tm': p['right_primer']['tm'],
         'self_complementarity': p['right_primer']['self_complementarity'],
         "self_3prime_complementarity": p['right_primer']['self_3prime_complementarity'],
         'amplicon_size': p['amplicon_size'],
         'amplicon_size_abs': p['amplicon_size_abs'],
         'label': f'Right {i + 1}', 'y': (i + 1) * 10}
        for i, p in enumerate(primers)]

    primers_data = pd.DataFrame(primers_data)

    # Définir un offset y fixe pour les exons
    y_offset = (len(primers) + 1) * 10  # Espacement dynamique basé sur le nombre de primers
    exons_data['y'] = [y_offset for _ in range(len(exons_data))]  # Toutes les positions y sont égales

    # Exons chart avec y fixe pour tous les exons et même hauteur
    exons_chart = alt.Chart(exons_data).mark_rect(stroke='darkblue', strokeWidth=2).encode(
        x=alt.X('start:Q', scale=alt.Scale(domain=[int(min(exons_data['start'])) - 1, int(max(exons_data['end']))]),
                title='Position (Exon + Intron; bp)' if normalization is False else "Position (Exon; bp)"),
        x2=alt.X2('end:Q'),
        y=alt.Y('y:Q', scale=alt.Scale(domain=[0, (len(primers) + 1) * 10])),  # Position y fixe
        color=alt.value('lightblue'),
        tooltip=['exon_number', 'start', 'end']
    ).interactive().properties(
        width=1500,
        height=600,
    )

    # Primers chart, maintenant colorié par jeu de primers
    primers_chart = alt.Chart(primers_data).mark_bar(height=10).encode(
        x=alt.X('position_start:Q' if normalization is True else 'position_start_abs:Q'),
        x2=alt.X2('position_end:Q' if normalization is True else 'position_end_abs:Q'),
        y=alt.Y('y:Q', title=None),
        color=alt.Color('pair:N', scale=alt.Scale(domain=[f'Pair {i + 1}' for i in range(len(primers))]),
                        title='Primers pairs'),
        tooltip=['pair', 'label', 'sequence', 'tm', 'length', 'position', 'self_complementarity', "self_3prime_complementarity", 'amplicon_size', 'position_abs', 'amplicon_size_abs'],
    )

    annotations_left = alt.Chart(primers_data).mark_text(
        align='left', dx=-40, fontSize=10, color='black'
    ).encode(
        x='position_start:Q' if normalization is True else 'position_start_abs:Q',
        y='y:Q',
        text='label:N'
    ).transform_filter(
        alt.datum.type == 'Left Primer'  # Filtrer pour les Left Primers
    )

    # Annotations pour les primers Right
    annotations_right = alt.Chart(primers_data).mark_text(
        align='right', dx=50, fontSize=10, color='black'
    ).encode(
        x='position_end:Q' if normalization is True else 'position_end_abs:Q',
        y='y:Q',
        text='label:N'
    ).transform_filter(
        alt.datum.type == 'Right Primer'  # Filtrer pour les Right Primers
    )

    # Combiner Exons et Primers
    final_chart = exons_chart + primers_chart + annotations_left + annotations_right

    final_chart.save("Test.html")
    return final_chart


# Page config
page_config()

st.subheader(':blue[Step 1] Promoter and Terminator Extractor')
colextract1, colextract2, colextract3 = st.columns([0.5, 1.5, 1.5], gap="small")

# Extraction of DNA sequence
with colextract1:
    st.info("💡 If you have a FASTA sequence, go to :blue[**Step 2**]")

    all_variants = {}
    upstream_entry = []

    # Gene ID
    st.markdown("🔹 :blue[**Step 1.1**] Gene ID:", help='NCBI gene name and NCBI gene ID allowed')
    gene_id_entry = st.text_area("🔹 :blue[**Step 1.1**] Gene ID:", value="PRKN\n351\nNM_003130.4",
                                 label_visibility='collapsed')
    gene_ids = gene_id_entry.strip().split("\n")

with colextract2:
    tab1, tab2 = st.tabs(['Default', 'Advance'])

    with tab1:
        # Species
        st.markdown("🔹 :blue[**Step 1.2**] Species of gene names and sliced variants:")
        col1, col2 = st.columns(2)
        # gr = col1.selectbox("Genome:", ["Current", "Previous"], index=0,
        #                     help='Example for Homo sapiens:\n\n"Current" is GRCh38\n\n"Previous" is GRCh37')
        species = col1.selectbox("Species:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)
        col2.markdown("")
        col2.markdown("")
        all_slice_form = col2.toggle(label='All variants')

        # Run Promoter Finder
        if st.button(f"🧬 :blue[**Step 1.5**] Extraction...", help='(~5sec/gene)'):
            response = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')

            ncbi_status = True if response.status_code == 200 else False

            if ncbi_status is True:
                with st.spinner('Please wait...'):
                    with colextract1:

                        pbar = st.progress(0,
                                           text='**:blue[Extract info...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                        for i, gene_id in enumerate(gene_ids):
                            pbar.progress(i / len(gene_ids),
                                          text=f'**:blue[Extract info... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                            all_variants_output, message = NCBIdna(gene_id, species, all_slice_forms=True if all_slice_form else False).find_sequences()
                            if "Error 200" not in all_variants_output:
                                pbar.progress((i + 1) / len(gene_ids),

                                              text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                                st.toast(f"**{gene_id}** from **{species}** extracted", icon='🧬')

                                all_variants.update(all_variants_output)

                            else:
                                st.error(message)
                                continue

                        st.session_state['all_variants'] = all_variants
                        st.success(f"Info extraction complete !")
                        st.toast(f"Info extraction complete !", icon='😊')

            elif ncbi_status is False:
                st.warning("⚠ NCBI servers are under maintenance or have an error")

    with tab2:
        # Advance mode extraction
        data_df = pd.DataFrame(
            {
                "Gene": gene_ids,
                "human": [False] * len(gene_ids),
                "mouse": [False] * len(gene_ids),
                "rat": [False] * len(gene_ids),
                "drosophila": [False] * len(gene_ids),
                "zebrafish": [False] * len(gene_ids)
            }
        )

        species_list = ['human', 'mouse', 'rat', 'drosophila', 'zebrafish']
        search_types = ['promoter', 'terminator']

        st.markdown('**🔹 :blue[Step 1.2]** Select species for all genes:',
                    help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')

        species1, species2, species3, species4, species5 = st.columns(5)

        with species1:
            all_human = st.toggle("Human")
        with species2:
            all_mouse = st.toggle("Mouse")
        with species3:
            all_rat = st.toggle("Rat")
        with species4:
            all_droso = st.toggle("Drosophila")
        with species5:
            all_zebra = st.toggle("Zebrafish")

        if all_human:
            data_df["human"] = True
        if all_mouse:
            data_df["mouse"] = True
        if all_rat:
            data_df["rat"] = True
        if all_droso:
            data_df["drosophila"] = True
        if all_zebra:
            data_df["zebrafish"] = True

        st.markdown('**🔹 :blue[Step 1.2]** On demand genes table',
                    help="Check the boxes for which you want to extract a sequence. Pay attention that the gene name is equivalent for each species. The choice of species is not available for gene IDs. Parameterize the table last, if you check the boxes above, it resets the whole table.")

        data_dff = st.data_editor(
            data_df,
            column_config={
                "human": st.column_config.CheckboxColumn(
                    "Human",
                    default=False,
                ),
                "mouse": st.column_config.CheckboxColumn(
                    "Mouse",
                    default=False,
                ),
                "rat": st.column_config.CheckboxColumn(
                    "Rat",
                    default=False,
                ),
                "drosophila": st.column_config.CheckboxColumn(
                    "Drosophila",
                    default=False,
                ),
                "zebrafish": st.column_config.CheckboxColumn(
                    "Zebrafish",
                    default=False,
                )
            },
            disabled=["Gene"],
            hide_index=True,
        )

        if st.button("🧬 :blue[**Step 1.4**] Extract sequences", help="(~5sec/seq)", key='Advance'):
            response = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')

            ncbi_status = True if response.status_code == 200 else False

            if ncbi_status is True:
                with colextract1:
                    pbar = st.progress(0,
                                       text='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                    for i, gene_info in enumerate(data_dff.itertuples(index=False)):
                        gene_id = gene_info.Gene
                        if gene_id.isdigit() or gene_id.startswith('XM_') or gene_id.startswith(
                                'NM_') or gene_id.startswith('XR_') or gene_id.startswith('NR_'):
                            for search_type in search_types:
                                if getattr(gene_info, f'{search_type}'):
                                    pbar.progress((i + 1) / len(data_dff),
                                                  text=f'**:blue[Extract sequence...**{gene_id}** from **{species}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                    all_variants_output, message = NCBIdna(gene_id).find_sequences()

                                    if "Error 200" not in all_variants_output:
                                        st.toast(f"**{gene_id}** from **{species}** extracted",
                                                 icon='🧬')

                                        all_variants.update(all_variants_output)
                                    else:
                                        st.error(message)
                                        continue

                        else:
                            for species in species_list:
                                for search_type in search_types:
                                    if getattr(gene_info, f'{species}') and getattr(gene_info,
                                                                                    f'{search_type}'):

                                        pbar.progress((i + 1) / len(data_dff),
                                                      text=f'**:blue[Extract sequence... **{gene_id}** from **{species.capitalize()}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                        all_variants_output, message = NCBIdna(gene_id, species).find_sequences()

                                        if "Error 200" not in all_variants_output:
                                            st.toast(f"**{gene_id}** from **{species}** extracted",
                                                     icon='🧬')
                                            all_variants.update(all_variants_output)
                                        else:
                                            st.error(message)
                                            continue

                    st.session_state['all_variants'] = all_variants
                    st.success(f"Info extraction complete !")
                    st.toast(f"Info extraction complete !", icon='😊')

            elif ncbi_status is False:
                st.warning("⚠ NCBI servers are under maintenance or have an error")

with colextract3:
    st.markdown("🔹 :blue[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')

    if 'all_variants' not in st.session_state:
        all_variants = ''
        st.session_state['all_variants'] = all_variants
    if len(st.session_state['all_variants']) > 0:
        output = st.json(st.session_state['all_variants'] if len(st.session_state['all_variants']) > 0 else {'Empty'}, expanded=False)

        current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        st.download_button(label="💾 Download (.json)", data=json.dumps(st.session_state['all_variants'], indent=4),
                           file_name=f"Sequences_{current_date_time}.json", mime="application/json")
    else:
        st.warning('You have not extracted any information')

# Zone de texte pour saisir une séquence FASTA
default_fasta = """>ExampleGene1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTA
>ExampleGene2
CGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""


fasta_input = st.text_area(
    "Enter a DNA sequence in FASTA format (required):",
    height=200,
    value=default_fasta
)


def parse_fasta_to_json(fasta_text):
    lines = fasta_text.strip().split("\n")
    variants = {}
    variant_id = 1
    gene_name = None
    sequence = ""

    for line in lines:
        if line.startswith(">"):
            if gene_name and sequence:
                variants[variant_id] = {
                    "gene_name": gene_name,
                    "sequence": sequence.upper(),
                    "normalized_exon_coords": [[0, len(sequence)]],
                }
                variant_id += 1
            gene_name = line[1:].strip()
            sequence = ""
        else:
            sequence += line.strip()

    if gene_name and sequence:
        variants[variant_id] = {
            "gene_name": gene_name,
            "sequence": sequence.upper(),
            "normalized_exon_coords": [[0, len(sequence)]],
        }

    return variants


if st.button("Add Sequence"):
    if fasta_input.strip():
        try:
            variants = parse_fasta_to_json(fasta_input)

            st.session_state['all_variants'] = variants
            st.success(f"{len(variants)} sequence(s) added successfully!")
            st.rerun()
        except Exception as e:
            st.error(f"Error parsing FASTA: {e}")
    else:
        st.error("Please provide a valid FASTA sequence.")

nb_primers = st.number_input('Number of primers', value=10, min_value=1, step=1)

if st.button('Run design primers'):
    try:
        primers_result = []

        if len(st.session_state['all_variants']) > 0:
            progress_bar = stqdm(enumerate(st.session_state['all_variants'].items()),
                                 total=len(st.session_state['all_variants']),
                                 desc="Initializing")

            for i, (variant, data) in progress_bar:
                # Mettre à jour la description à chaque itération
                progress_bar.set_description(f"Designing primers for {variant} - {data['gene_name']}")

                primers = NCBIdna.design_primers(variant, data['gene_name'], data['sequence'], data['normalized_exon_coords'], nb_primers)

                if len(primers) > 0:
                    for idx, primer_set in enumerate(primers):
                        primers_result.append({
                            'Gene': str(variant) + data['gene_name'],
                            'Pair': idx + 1,
                            'Product Size (bp)': primer_set['amplicon_size'],
                            "For. Pr.(5'->3')": primer_set['left_primer']['sequence'],
                            'For. Len. (bp)': primer_set['left_primer']['length'],
                            'For. Pos.': primer_set['left_primer']['position'],
                            'For. Pos. Abs. (bp)': primer_set['left_primer']['position_abs'],
                            'For. Tm (°C)': primer_set['left_primer']['tm'],
                            'For. GC%': primer_set['left_primer']['gc_percent'],
                            'For. Self Compl.': primer_set['left_primer']['self_complementarity'],
                            "For. Self 3' Compl.": primer_set['left_primer']['self_3prime_complementarity'],
                            "Rev. Pr.(5'->3')": primer_set['right_primer']['sequence'],
                            'Rev. Len. (bp)': primer_set['right_primer']['length'],
                            'Rev. Pos.': primer_set['right_primer']['position'],
                            'Rev. Pos. Abs.': primer_set['right_primer']['position_abs'],
                            'Rev. Tm (°C)': primer_set['right_primer']['tm'],
                            'Rev. GC%': primer_set['right_primer']['gc_percent'],
                            'Rev. Self Compl.': primer_set['right_primer']['self_complementarity'],
                            "Rev. Self 3' Compl.": primer_set['right_primer']['self_3prime_complementarity'],
                            'Product Size Abs. (bp)': primer_set['amplicon_size_abs']
                        })

                    with st.expander(f'Primers graph location for {variant} {data["gene_name"]}', expanded=False):
                        if len(data['normalized_exon_coords']) > 1:
                            st.altair_chart(graphique(data['normalized_exon_coords'], primers), theme=None,
                                            use_container_width=True,
                                            key=f"{variant}_exon_intron")
                        st.altair_chart(graphique(data['normalized_exon_coords'], primers, True), theme=None,
                                        use_container_width=True, key=f"{variant}_exon")

                    st.toast(f"Primers designed for {variant} {data['gene_name']}!")
                else:
                    st.warning(f"No primers were designed for {variant} {data['gene_name']}")

            if primers_result:
                st.dataframe(pd.DataFrame(primers_result), hide_index=True)

                csv_file = pd.DataFrame(primers_result).to_csv(index=False)
                excel_file = io.BytesIO()
                pd.DataFrame(primers_result).to_excel(excel_file, index=False, sheet_name='Sheet1')
                excel_file.seek(0)

                download_button1, download_button2 = st.columns(2, gap='small')
                current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                download_button1.download_button("💾 Download table (.xlsx)", excel_file,
                                   file_name=f'LabmasterDP_{current_date_time}.xlsx',
                                   mime="application/vnd.ms-excel", key='download-excel')
                download_button2.download_button(label="💾 Download table (.csv)", data=csv_file,
                                   file_name=f"LabmasterDP_{current_date_time}.csv", mime="text/csv")
    except Exception as e:
        print(e)



ucsc_species = {
    "Homo sapiens": {"uscs_species": "Human"},
    "Mus musculus": {"uscs_species": "Mouse"},
    "Anopheles gambiae": {"uscs_species": "A. gambiae"},
    "Apis mellifera": {"uscs_species": "A. mellifera"},
    "Xenopus laevis": {"uscs_species": "African clawed frog"},
    "Vicugna pacos": {"uscs_species": "Alpaca"},
    "Alligator mississippiensis": {"uscs_species": "American alligator"},
    "Dasypus novemcinctus": {"uscs_species": "Armadillo"},
    "Gadus morhua": {"uscs_species": "Atlantic cod"},
    "Papio anubis": {"uscs_species": "Baboon"},
    "Bison bison": {"uscs_species": "Bison"},
    "Pan paniscus": {"uscs_species": "Bonobo"},
    "Apteryx australis": {"uscs_species": "Brown kiwi"},
    "Melopsittacus undulatus": {"uscs_species": "Budgerigar"},
    "Otolemur garnettii": {"uscs_species": "Bushbaby"},
    "Caenorhabditis brenneri": {"uscs_species": "C. brenneri"},
    "Caenorhabditis briggsae": {"uscs_species": "C. briggsae"},
    "Caenorhabditis elegans": {"uscs_species": "C. elegans"},
    "Ciona intestinalis": {"uscs_species": "C. intestinalis"},
    "Caenorhabditis japonica": {"uscs_species": "C. japonica"},
    "Caenorhabditis remanei": {"uscs_species": "C. remanei"},
    "Felis catus": {"uscs_species": "Cat"},
    "Gallus gallus": {"uscs_species": "Chicken"},
    "Pan troglodytes": {"uscs_species": "Chimp"},
    "Cricetulus griseus": {"uscs_species": "Chinese hamster"},
    "Manis pentadactyla": {"uscs_species": "Chinese pangolin"},
    "Latimeria chalumnae": {"uscs_species": "Coelacanth"},
    "Bos taurus": {"uscs_species": "Cow"},
    "Macaca fascicularis": {"uscs_species": "Crab-eating macaque"},
    "Drosophila ananassae": {"uscs_species": "D. ananassae"},
    "Drosophila erecta": {"uscs_species": "D. erecta"},
    "Drosophila grimshawi": {"uscs_species": "D. grimshawi"},
    "Drosophila melanogaster": {"uscs_species": "D. melanogaster"},
    "Drosophila mojavensis": {"uscs_species": "D. mojavensis"},
    "Drosophila persimilis": {"uscs_species": "D. persimilis"},
    "Drosophila pseudoobscura": {"uscs_species": "D. pseudoobscura"},
    "Drosophila sechellia": {"uscs_species": "D. sechellia"},
    "Drosophila simulans": {"uscs_species": "D. simulans"},
    "Drosophila virilis": {"uscs_species": "D. virilis"},
    "Drosophila yakuba": {"uscs_species": "D. yakuba"},
    "Canis lupus familiaris": {"uscs_species": "Dog"},
    "Delphinidae": {"uscs_species": "Dolphin"},
    "Zaire ebolavirus": {"uscs_species": "Ebola virus"},
    "Loxodonta africana": {"uscs_species": "Elephant"},
    "Callorhinchus milii": {"uscs_species": "Elephant shark"},
    "Mustela putorius furo": {"uscs_species": "Ferret"},
    "Takifugu rubripes": {"uscs_species": "Fugu"},
    "Thamnophis sirtalis": {"uscs_species": "Garter snake"},
    "Hylobates": {"uscs_species": "Gibbon"},
    "Aquila chrysaetos": {"uscs_species": "Golden eagle"},
    "Rhinopithecus roxellana": {"uscs_species": "Golden snub-nosed monkey"},
    "Gorilla gorilla": {"uscs_species": "Gorilla"},
    "Chlorocebus sabaeus": {"uscs_species": "Green monkey"},
    "Cavia porcellus": {"uscs_species": "Guinea pig"},
    "Neomonachus schauinslandi": {"uscs_species": "Hawaiian monk seal"},
    "Erinaceus europaeus": {"uscs_species": "Hedgehog"},
    "Equus caballus": {"uscs_species": "Horse"},
    "Dipodomys": {"uscs_species": "Kangaroo rat"},
    "Petromyzon marinus": {"uscs_species": "Lamprey"},
    "Branchiostoma": {"uscs_species": "Lancelet"},
    "Sceloporus": {"uscs_species": "Lizard"},
    "Cynocephalus variegatus": {"uscs_species": "Malayan flying lemur"},
    "Trichechus manatus": {"uscs_species": "Manatee"},
    "Callithrix jacchus": {"uscs_species": "Marmoset"},
    "Oryzias latipes": {"uscs_species": "Medaka"},
    "Geospiza fortis": {"uscs_species": "Medium ground finch"},
    "Pteropus": {"uscs_species": "Megabat"},
    "Myotis lucifugus": {"uscs_species": "Little brown bat"},
    "Balaenoptera acutorostrata": {"uscs_species": "Minke whale"},
    "Microcebus murinus": {"uscs_species": "Mouse lemur"},
    "Heterocephalus glaber": {"uscs_species": "Naked mole-rat"},
    "Oreochromis niloticus": {"uscs_species": "Nile tilapia"},
    "Monodelphis domestica": {"uscs_species": "Opossum"},
    "Pongo pygmaeus": {"uscs_species": "Orangutan"},
    "Pristionchus pacificus": {"uscs_species": "P. pacificus"},
    "Chrysemys picta": {"uscs_species": "Painted turtle"},
    "Ailuropoda melanoleuca": {"uscs_species": "Panda"},
    "Sus scrofa": {"uscs_species": "Pig"},
    "Ochotona princeps": {"uscs_species": "Pika"},
    "Ornithorhynchus anatinus": {"uscs_species": "Platypus"},
    "Nasalis larvatus": {"uscs_species": "Proboscis monkey"},
    "Oryctolagus cuniculus": {"uscs_species": "Rabbit"},
    "Rattus norvegicus": {"uscs_species": "Rat"},
    "Macaca mulatta": {"uscs_species": "Rhesus"},
    "Procavia capensis": {"uscs_species": "Rock hyrax"},
    "Saccharomyces cerevisiae": {"uscs_species": "S. cerevisiae"},
    "Strongylocentrotus purpuratus": {"uscs_species": "S. purpuratus"},
    "Aplysia californica": {"uscs_species": "Sea hare"},
    "Enhydra lutris nereis": {"uscs_species": "Southern sea otter"},
    "Ovis aries": {"uscs_species": "Sheep"},
    "Sorex araneus": {"uscs_species": "Shrew"},
    "Bradypus": {"uscs_species": "Sloth"},
    "Sciurus": {"uscs_species": "Squirrel"},
    "Saimiri sciureus": {"uscs_species": "Squirrel monkey"},
    "Gasterosteus aculeatus": {"uscs_species": "Stickleback"},
    "Tarsius syrichta": {"uscs_species": "Tarsier"},
    "Sarcophilus harrisii": {"uscs_species": "Tasmanian devil"},
    "Tenrec ecaudatus": {"uscs_species": "Tenrec"},
    "Tetraodon nigroviridis": {"uscs_species": "Tetraodon"},
    "Nanorana parkeri": {"uscs_species": "Tibetan frog"},
    "Tupaia": {"uscs_species": "Tree shrew"},
    "Meleagris gallopavo": {"uscs_species": "Turkey"},
    "Severe acute respiratory syndrome coronavirus 2": {"uscs_species": "SARS-CoV-2"},
    "Monkeypox virus": {"uscs_species": "Monkeypox virus"},
    "Macropus": {"uscs_species": "Wallaby"},
    "Ceratotherium simum": {"uscs_species": "White rhinoceros"},
    "Xenopus tropicalis": {"uscs_species": "X. tropicalis"},
    "Taeniopygia guttata": {"uscs_species": "Zebra finch"},
    "Danio rerio": {"uscs_species": "Zebrafish"}
}


