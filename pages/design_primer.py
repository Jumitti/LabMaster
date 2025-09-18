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


ucsc_species = [
    "None",
    "Homo sapiens",
    "Mus musculus",
    "Anopheles gambiae",
    "Apis mellifera",
    "Xenopus laevis",
    "Vicugna pacos",
    "Alligator mississippiensis",
    "Dasypus novemcinctus",
    "Gadus morhua",
    "Papio anubis",
    "Bison bison",
    "Pan paniscus",
    "Apteryx australis",
    "Melopsittacus undulatus",
    "Otolemur garnettii",
    "Caenorhabditis brenneri",
    "Caenorhabditis briggsae",
    "Caenorhabditis elegans",
    "Ciona intestinalis",
    "Caenorhabditis japonica",
    "Caenorhabditis remanei",
    "Felis catus",
    "Gallus gallus",
    "Pan troglodytes",
    "Cricetulus griseus",
    "Manis pentadactyla",
    "Latimeria chalumnae",
    "Bos taurus",
    "Macaca fascicularis",
    "Drosophila ananassae",
    "Drosophila erecta",
    "Drosophila grimshawi",
    "Drosophila melanogaster",
    "Drosophila mojavensis",
    "Drosophila persimilis",
    "Drosophila pseudoobscura",
    "Drosophila sechellia",
    "Drosophila simulans",
    "Drosophila virilis",
    "Drosophila yakuba",
    "Canis lupus familiaris",
    "Delphinidae",
    "Zaire ebolavirus",
    "Loxodonta africana",
    "Callorhinchus milii",
    "Mustela putorius furo",
    "Takifugu rubripes",
    "Thamnophis sirtalis",
    "Hylobates",
    "Aquila chrysaetos",
    "Rhinopithecus roxellana",
    "Gorilla gorilla",
    "Chlorocebus sabaeus",
    "Cavia porcellus",
    "Neomonachus schauinslandi",
    "Erinaceus europaeus",
    "Equus caballus",
    "Dipodomys",
    "Petromyzon marinus",
    "Branchiostoma",
    "Sceloporus",
    "Cynocephalus variegatus",
    "Trichechus manatus",
    "Callithrix jacchus",
    "Oryzias latipes",
    "Geospiza fortis",
    "Pteropus",
    "Myotis lucifugus",
    "Balaenoptera acutorostrata",
    "Microcebus murinus",
    "Heterocephalus glaber",
    "Oreochromis niloticus",
    "Monodelphis domestica",
    "Pongo pygmaeus",
    "Pristionchus pacificus",
    "Chrysemys picta",
    "Ailuropoda melanoleuca",
    "Sus scrofa",
    "Ochotona princeps",
    "Ornithorhynchus anatinus",
    "Nasalis larvatus",
    "Oryctolagus cuniculus",
    "Rattus norvegicus",
    "Macaca mulatta",
    "Procavia capensis",
    "Saccharomyces cerevisiae",
    "Strongylocentrotus purpuratus",
    "Aplysia californica",
    "Enhydra lutris nereis",
    "Ovis aries",
    "Sorex araneus",
    "Bradypus",
    "Sciurus",
    "Saimiri sciureus",
    "Gasterosteus aculeatus",
    "Tarsius syrichta",
    "Sarcophilus harrisii",
    "Tenrec ecaudatus",
    "Tetraodon nigroviridis",
    "Nanorana parkeri",
    "Tupaia",
    "Meleagris gallopavo",
    "Severe acute respiratory syndrome coronavirus 2",
    "Monkeypox virus",
    "Macropus",
    "Ceratotherium simum",
    "Xenopus tropicalis",
    "Taeniopygia guttata",
    "Danio rerio"
]


def parse_fasta_to_json(fasta_text):
    lines = fasta_text.strip().split("\n")
    variants = {}
    variant_id = 1
    sequence = ""

    # Dictionnaire en minuscule pour une recherche insensible à la casse
    ucsc_species_lower = {key.lower(): key for key in ucsc_species}

    for line in lines:
        if line.startswith(">"):
            # Sauvegarder l'entrée précédente si elle existe
            if sequence:
                variants[variant_id]["sequence"] = sequence.upper()
                variants[variant_id]["normalized_exon_coords"] = [[0, len(sequence)]]
                variant_id += 1

            # Extraire le nom du gène et l'espèce s'il y a un "|"
            parts = line[1:].strip().split("|")
            gene_name = parts[0].strip()
            species_name = parts[1].strip() if len(parts) > 1 else ""

            # Vérifier si l'espèce correspond à une clé UCSC
            species = ucsc_species_lower.get(species_name.lower(), None)

            # Créer une nouvelle entrée
            variants[variant_id] = {
                "gene_name": gene_name,
                "species": species,
                "sequence": "",
                "normalized_exon_coords": [],
            }

            sequence = ""  # Réinitialiser la séquence pour la prochaine entrée

        else:
            sequence += line.strip()

    # Ajouter la dernière entrée après la boucle
    if sequence:
        variants[variant_id]["sequence"] = sequence.upper()
        variants[variant_id]["normalized_exon_coords"] = [[0, len(sequence)]]

    return variants


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

st.subheader('Genes info extraction')
colextract1, colextract2, colextract3 = st.columns([0.5, 1.5, 1.5], gap="small")

# Extraction of DNA sequence
with colextract1:

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
        if st.button(f"🧬 :blue[**Step 1.3**] Extraction...", help='(~5sec/gene)'):
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
    default_fasta = """>ExampleGene1 | Homo sapiens\ntacgcaatgtcatgactgcgtttatagatagataaaaagcgtgcgattactaaacgcggatggcgtgcgcactatttcatcggttctgaaatctccatccaatcaaccttactcagacagctcccccgtgacacgggctaccacattcaggtggcttgtaatacatgggtataacatcaatagttcgtgccgcaatacttcgcgggggtacgggtaagtgacgaaagaagtaactctcccactcggagaatctacggtagttgcgcgtttttaattttcatctttgtcctgccagcaatgtacacaccgcaaagtctgtccaagtgcatgctagaccgggtgtgcaccctagggtagagcacggagttgatttcgggcgtgagatcaaggccaaggaggaagtaagcatcgtatctctgtctaatcattgcaggaagggtgcacagcttaggttccccaacaatgcttttagcatgatagctgtctctttgtggactgta\n>ExampleGene2\nttcctctttcaccagctttccatccccgcgacatggcggatcaaaactctggcaaagattaccagtcgaaggcatctcgagatggagatggtaagtttttgtcatacgacccaaacccggaggagacacgttagaaaatccacgacttcttcgaagactaagtggatgagtacaggtcgggaagagtcgactaccctaggatcccgcgtgcggtctacatgtcatgatcctccatgggcccaggccccgtagtgcgactgcggttaattgcatctacgaattttacacttgcgtttaagaccggacgccgggtttctaagtaaaagtttggctatcgacattatttttttaggggcaccgtatcggattccaatgggtggctggattctcagtgaatctccgtagttcgggaaatcactcaggaatgctaatcatccagaatggaaacgtggtaaaagactgccctgcttccctcttttacctcaagaaacaggggcggg"""

    st.markdown('**🔹 :blue[Alternative 1.1]** Enter a DNA sequence in FASTA format (required):',)
    with st.expander("📌 How to write a FASTA file:"):
        st.markdown("Each entry should follow this structure:")

        st.code(""">GeneName or ID or other | Species\nATGCATGCATGC""", language="plaintext")

        st.markdown("For example:")
        st.code(""">PRKN | Homo sapiens\nATGCGTGCATGCGTGCATGCGTGC""", language="plaintext")

    fasta_input = st.text_area(
        "Enter a DNA sequence in FASTA format (required):",
        height=200,
        value=default_fasta, label_visibility="collapsed"
    )

    if st.button("**🔹 :blue[Alternative 1.2]** Add your sequences"):
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

st.markdown("🔹 :blue[**Sequences added**]")
st.markdown("You can change the species if necessary for the validation of primers via UCSC PCR in-Silico. Use 'None' if species isn't known")


if 'all_variants' not in st.session_state:
    st.session_state['all_variants'] = {}

if st.session_state['all_variants']:
    df = pd.DataFrame.from_dict(st.session_state['all_variants'], orient="index").reset_index()
    df.rename(columns={"index": "N°/Name"}, inplace=True)

    if "species" not in df.columns:
        df["species"] = ""

    column_order = list(df.columns)
    if "species" in column_order and "exon_coords" in column_order:
        column_order.remove("species")
        exon_index = column_order.index("exon_coords")
        column_order.insert(exon_index, "species")
        df = df[column_order]

    column_config = {
        "species": st.column_config.SelectboxColumn("Species", options=ucsc_species, required=True)
    }
    for col in df.columns:
        if col != "species":
            column_config[col] = st.column_config.Column(col, disabled=True)

    df = st.data_editor(df, column_config=column_config, use_container_width=True, hide_index=True)

    updated_json = df.set_index("N°/Name").to_dict(orient="index")
    st.json(updated_json, expanded=False)
    current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    st.download_button(
        label="💾 Download (.json)",
        data=json.dumps(updated_json, indent=4),
        file_name=f"Sequences_{current_date_time}.json",
        mime="application/json"
    )
else:
    st.warning('You have not extracted any information')

if "min_amplicon_size" not in st.session_state:
    st.session_state["min_amplicon_size"] = 80
if "max_amplicon_size" not in st.session_state:
    st.session_state["max_amplicon_size"] = 250


nb_primers = st.number_input('Number of primers', value=10, min_value=1, step=1)
st.session_state["min_amplicon_size"] = st.number_input('Minimum amplicon size', value=60, min_value=10,
                                                        max_value=st.session_state["max_amplicon_size"] - 10, step=1)

st.session_state["max_amplicon_size"] = st.number_input('Maximum amplicon size', value=250,
                                                        min_value=st.session_state["min_amplicon_size"] + 10, step=1)

if 'ucsc_validation' not in st.session_state:
    st.session_state["ucsc_validation"] = False
if 'only_validated' not in st.session_state:
    st.session_state["only_validated"] = False

st.session_state["ucsc_validation"] = st.toggle("UCSC validation", False)
if st.session_state["ucsc_validation"] is True:
    st.session_state["only_validated"] = st.toggle("Display only validated", False)


if st.button('Run design primers'):
    try:
        primers_result = []

        if len(st.session_state['all_variants']) > 0:
            progress_bar = stqdm(enumerate(st.session_state['all_variants'].items()),
                                 total=len(st.session_state['all_variants']),
                                 desc="Initializing")
            for i, (variant, data) in progress_bar:
                progress_bar.set_description(f"Designing primers for {variant} - {data['gene_name']}")

                primers = NCBIdna.design_primers(variant, data['gene_name'], data['species'], data['sequence'], data['normalized_exon_coords'], nb_primers, [st.session_state["min_amplicon_size"], st.session_state["max_amplicon_size"]], st.session_state["ucsc_validation"], st.session_state["only_validated"])

                if len(primers) > 0:
                    for idx, primer_set in enumerate(primers):
                        primers_result.append({
                            'Gene': str(variant) + data['gene_name'],
                            'Pair': idx + 1,
                            'Product Size (bp)': primer_set['amplicon_size'],
                            'Validated': primer_set['validation_relative'],
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
                            'Product Size Abs. (bp)': primer_set['amplicon_size_abs'],
                            'Validated Abs.': primer_set['validation_absolute'],
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
                all_columns = list(primers_result[0].keys())
                abs_columns = [col for col in all_columns if "Abs." in col]
                other_columns = [col for col in all_columns if col not in abs_columns]
                ordered_columns = other_columns + abs_columns
                df = pd.DataFrame(primers_result)[ordered_columns]
                st.dataframe(df, hide_index=True)

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


