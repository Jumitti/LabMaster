import requests
import xml.etree.ElementTree as ET
import time
import random
import datetime
import json

from Bio.Seq import Seq
from Bio import Entrez
import primer3
from pages.tfinder import NCBIdna
import altair as alt
import pandas as pd
from tqdm import tqdm
from utils.page_config import page_config
import streamlit as st


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
                     'position': p['left_primer']['position'],
                     'position_start': p['left_primer']['position'][0],
                     'position_end': p['left_primer']['position'][1],
                     'position_abs': p['left_primer']['position_abs'],
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
         'position': p['right_primer']['position'],
         'position_start': p['right_primer']['position'][0],
         'position_end': p['right_primer']['position'][1],
         'position_abs': p['right_primer']['position_abs'],
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

    max_end = max(exons_data['end'])  # Cette ligne récupère la valeur maximale de la position "end"

    # Définir un offset y fixe pour les exons
    y_offset = (len(primers) + 1) * 10  # Espacement dynamique basé sur le nombre de primers
    exons_data['y'] = [y_offset for _ in range(len(exons_data))]  # Toutes les positions y sont égales

    # Exons chart avec y fixe pour tous les exons et même hauteur
    exons_chart = alt.Chart(exons_data).mark_rect(stroke='darkblue', strokeWidth=2).encode(
        x=alt.X('start:Q',
                title='Position (Exon + Intron; bp)' if normalization is False else "Position (Exon; bp)"),
        x2=alt.X2('end:Q'),
        y=alt.Y('y:Q', scale=alt.Scale(domain=[0, (len(primers) + 1) * 10])),  # Position y fixe
        color=alt.value('lightblue'),
        tooltip=['exon_number', 'start', 'end']
    ).interactive().properties(
        width=1920,
        height=1080,
    )

    # Primers chart, maintenant colorié par jeu de primers
    primers_chart = alt.Chart(primers_data).mark_bar(height=10).encode(
        x=alt.X('position_start:Q' if normalization is True else 'position_start_abs:Q'),
        x2=alt.X2('position_end:Q' if normalization is True else 'position_end_abs:Q'),
        y=alt.Y('y:Q', title=None),  # Utilise la colonne 'y' pour espacer les points
        color=alt.Color('pair:N', scale=alt.Scale(domain=[f'Pair {i + 1}' for i in range(len(primers))]),
                        title='Primers pairs'),
        tooltip=['pair', 'label', 'sequence', 'tm', 'position', 'position_abs', 'tm', 'self_complementarity', "self_3prime_complementarity", 'amplicon_size', 'amplicon_size_abs']
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

    # Sauvegarder ou afficher
    final_chart.save(f'sequence_visualization{"_adjusted" if normalization is True else ""}.html')
    print(f"Graphique sauvegardé sous 'sequence_visualization{'_adjusted' if normalization is True else ''}.html'. Ouvrez ce fichier pour le visualiser.")


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
    output = st.json(st.session_state['all_variants'], expanded=False)

    current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    st.download_button(label="💾 Download (.json)", data=json.dumps(st.session_state['all_variants'], indent=4),
                       file_name=f"Sequences_{current_date_time}.json", mime="application/json")

'''
# Exemple d'utilisation
entrez_id = "4843"  # ID du gène BRCA1
# Récupérer toutes les variantes et leurs coordonnées d'exons
all_variants, message = NCBIdna.all_variant(entrez_id)
print("All variants:", all_variants)

# Exemple d'extraction de séquence à partir des coordonnées des exons
for variant, gene_name, chromosome, exon_coords, normalized_coords, species_API in all_variants:
    if exon_coords:
        # Utiliser le premier start et le dernier end pour extraire la séquence
        first_start = exon_coords[0][0]
        last_end = exon_coords[-1][1]

        # Récupérer la séquence entre ces coordonnées
        sequence = NCBIdna.get_dna_sequence(gene_name, chromosome, first_start, last_end)
        print(f"Sequence extracted for {gene_name}: {sequence}")

        print(normalized_coords)


        def design_primers(sequence, exons, nb_primers):
            simplified_sequence = "".join(sequence[start:end] for start, end in exons)

            primer3_input = {
                'SEQUENCE_ID': 'primer_in_exons',
                'SEQUENCE_TEMPLATE': simplified_sequence,
            }

            primer3_params = {
                'PRIMER_TASK': 'generic',  # Generic primer generation

                # Settings for primers (size and content)
                'PRIMER_OPT_SIZE': 20,  # Optimal primer size
                'PRIMER_MIN_SIZE': 18,  # Minimum primer size
                'PRIMER_MAX_SIZE': 24,  # Maximum primer size
                'PRIMER_OPT_TM': 60.0,  # Optimal melting temperature (°C)
                'PRIMER_MIN_TM': 57.0,  # Minimum melting temperature (°C)
                'PRIMER_MAX_TM': 63.0,  # Maximum melting temperature (°C)
                'PRIMER_MIN_GC': 45.0,  # Minimum GC percentage (%)
                'PRIMER_MAX_GC': 55.0,  # Maximum GC percentage (%)
                'PRIMER_GC_CLAMP': 0,  # GC clamping at end 3' (minimum number of G/C)
                'PRIMER_MAX_POLY_X': 5,  # Maximum number of repeated bases (ex: AAAAA)

                # Parameters for stability at the 3' end
                'PRIMER_MAX_END_STABILITY': 9.0,  # Maximum end stability 3'

                # Secondary alignment (Thermodynamic model)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 70.0,  # Bad template match (primer pairs)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH_TMPL': 40.0,  # Bad match for single primer
                'PRIMER_MAX_SELF_ANY_TH': 45.0,  # Internal matching (all sites, thermodynamics)
                'PRIMER_MAX_SELF_END_TH': 35.0,  # Internal pairing (3' end, thermodynamic)
                'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,  # Primer pairing (all sites, thermodynamics)
                'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,  # Pairing between primers (3' end, thermodynamics)
                'PRIMER_MAX_HAIRPIN_TH': 24.0,  # Maximum free energy for hairpins

                # Secondary alignment (Old model)
                'PRIMER_MAX_TEMPLATE_MISPRIMING': 24.0,  # Bad template matching (primer pairs, classic)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TMPL': 12.0,  # Bad match for single primer (classic)
                'PRIMER_MAX_SELF_ANY': 8.0,  # Internal pairing (all sites, classic)
                'PRIMER_MAX_SELF_END': 3.0,  # Internal pairing (3' end, classic)
                'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,  # Primer pairing (all sites, classic)
                'PRIMER_PAIR_MAX_COMPL_END': 3.0,  # Pairing between primers (3' end, classic)

                # Search for secondary alignments
                'PRIMER_THERMODYNAMIC_ALIGNMENT': 1,  # Use thermodynamic model
                'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1, # Also align with thermodynamic model (maybe slow)

                # General settings for pairs
                'PRIMER_NUM_RETURN': 0,  # Maximum number of pairs returned
                'PRIMER_PRODUCT_SIZE_RANGE': [[80, 250]],  # Product size range

                # Enable selection of internal hybridization oligos
                'PRIMER_PICK_INTERNAL_OLIGO': 0,  # 1 to enable, 0 to disable

                # Size parameters for internal oligos
                'PRIMER_INTERNAL_MIN_SIZE': 18,  # Minimum size
                'PRIMER_INTERNAL_OPT_SIZE': 20,  # Optimal size
                'PRIMER_INTERNAL_MAX_SIZE': 27,  # Maximum size

                # Melting temperature (Tm) for internal oligos
                'PRIMER_INTERNAL_MIN_TM': 57.0,  # Minimum Tm (°C)
                'PRIMER_INTERNAL_OPT_TM': 60.0,  # Optimal Tm (°C)
                'PRIMER_INTERNAL_MAX_TM': 63.0,  # Maximum Tm (°C)

                # GC percentage for internal oligos
                'PRIMER_INTERNAL_MIN_GC': 20.0,  # Minimum GC percentage
                'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,  # Optimal GC percentage
                'PRIMER_INTERNAL_MAX_GC': 80.0,  # Maximum GC percentage

                # Concentration of monovalent cations (e.g.: Na+)
                'PRIMER_MONOVALENT_CATION_CONC': 50.0,  # In mM (default is 50 mM)

                # Concentration of divalent cations (e.g. Mg2+)
                'PRIMER_DIVALENT_CATION_CONC': 1.5,  # In mM (default is 1.5 mM)

                # Concentration of dNTPs
                'PRIMER_DNTP_CONC': 0.6,  # In mM (default is 0.6 mM)

                # Salt correction formula (using Santa Lucia 1998 formula)
                'PRIMER_SALT_CORRECTION': 1,  # 1 to use the Santa Lucia 1998 correction

                # Thermodynamic parameters (table of parameters to calculate Tm)
                'PRIMER_THERMODYNAMIC_PARAMETERS': 'SantaLucia1998',  # Use Santa Lucia 1998 table

                # Concentration of the oligonucleotide for the init
                'PRIMER_ANN_Oligo_CONC': 50.0,
            }

            primers = []

            # Calcul des longueurs cumulées des exons pour convertir en positions absolues
            exon_lengths = [end - start for start, end in exons]
            cumulative_lengths = [0] + list(cumsum(exon_lengths))

            # Combinaisons d'exons
            exon_pairs = [(i, j) for i in range(len(exons)) for j in range(i + 1, len(exons))]

            # Suivi des séquences déjà rencontrées pour éviter les doublons
            seen_primers = set()

            # Initialisez la barre de progression
            with tqdm(total=nb_primers, desc="Generating primers", unit="primer") as pbar:
                while len(primers) < nb_primers:
                    # Relance la recherche avec plus de résultats si nécessaire
                    primer3_params['PRIMER_NUM_RETURN'] += 1

                    for i, j in exon_pairs:
                        if len(primers) >= nb_primers:
                            break

                        # Définir les exons pour cette paire
                        exon1_start, exon1_end = exons[i]
                        exon2_start, exon2_end = exons[j]

                        # Simplifier les positions
                        simplified_start1 = cumulative_lengths[i]
                        simplified_end1 = simplified_start1 + (exon1_end - exon1_start)
                        simplified_start2 = cumulative_lengths[j]
                        simplified_end2 = simplified_start2 + (exon2_end - exon2_start)

                        # Taille du produit
                        product_size = simplified_end2 - simplified_start1

                        # Vérifier si la taille est valide
                        if 80 <= product_size <= 250:
                            primer3_input['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [
                                simplified_start1, simplified_end1 - simplified_start1,
                                simplified_start2, simplified_end2 - simplified_start2
                            ]

                            # Appeler Primer3 pour concevoir les amorces
                            primer_results = primer3.bindings.design_primers(primer3_input, primer3_params)

                            # Vérifier si des amorces ont été générées
                            if 'PRIMER_PAIR_NUM_RETURNED' in primer_results and primer_results[
                                'PRIMER_PAIR_NUM_RETURNED'] > 0:
                                for k in range(primer_results['PRIMER_PAIR_NUM_RETURNED']):
                                    if len(primers) >= nb_primers:
                                        break

                                    left_key = f'PRIMER_LEFT_{k}_SEQUENCE'
                                    right_key = f'PRIMER_RIGHT_{k}_SEQUENCE'

                                    if left_key in primer_results and right_key in primer_results:
                                        # Récupérer les séquences des amorces
                                        left_seq = primer_results.get(left_key, 'N/A')
                                        right_seq = primer_results.get(right_key, 'N/A')

                                        # Créer une clé unique pour détecter les doublons
                                        primer_key = (left_seq, right_seq)
                                        if primer_key in seen_primers:
                                            continue  # Passer si les amorces sont déjà enregistrées

                                        # Calculer les positions
                                        left_position = primer_results.get(f'PRIMER_LEFT_{k}')[0]
                                        right_position = primer_results.get(f'PRIMER_RIGHT_{k}')[0]

                                        left_absolute = convert_to_absolute(left_position, exons, cumulative_lengths)
                                        right_absolute = convert_to_absolute(right_position, exons, cumulative_lengths)

                                        amplicon_size = right_position - left_position + 1
                                        amplicon_size_abs = right_absolute - left_absolute + 1

                                        primers.append({
                                            'left_primer': {
                                                'sequence': left_seq,
                                                'length': len(left_seq),
                                                'position': (left_position, left_position + len(left_seq)),
                                                'position_abs': (left_absolute, left_absolute + len(left_seq)),
                                                'tm': primer_results.get(f'PRIMER_LEFT_{k}_TM', 'N/A'),
                                                'gc_percent': primer_results.get(f'PRIMER_LEFT_{k}_GC_PERCENT', 'N/A'),
                                                'self_complementarity': primer_results.get(f'PRIMER_LEFT_{k}_SELF_ANY_TH',
                                                                                           'N/A'),
                                                'self_3prime_complementarity': primer_results.get(
                                                    f'PRIMER_LEFT_{k}_SELF_END_TH', 'N/A'),
                                                'exon_junction': None
                                            },
                                            'right_primer': {
                                                'sequence': right_seq,
                                                'length': len(right_seq),
                                                'position': (right_position, right_position + len(right_seq)),
                                                'position_abs': (right_absolute, right_absolute + len(right_seq)),
                                                'tm': primer_results.get(f'PRIMER_RIGHT_{k}_TM', 'N/A'),
                                                'gc_percent': primer_results.get(f'PRIMER_RIGHT_{k}_GC_PERCENT', 'N/A'),
                                                'self_complementarity': primer_results.get(f'PRIMER_RIGHT_{k}_SELF_ANY_TH',
                                                                                           'N/A'),
                                                'self_3prime_complementarity': primer_results.get(
                                                    f'PRIMER_RIGHT_{k}_SELF_END_TH', 'N/A'),
                                                'template_strand': 'Minus',
                                                'exon_junction': None
                                            },
                                            'amplicon_size': amplicon_size,
                                            'amplicon_size_abs': amplicon_size_abs
                                        })

                                        pbar.update(1)
                                        seen_primers.add(primer_key)

            return primers


        def cumsum(iterable):
            total = 0
            for value in iterable:
                total += value
                yield total


        def convert_to_absolute(position, exons, cumulative_lengths):
            for i, (start, end) in enumerate(exons):
                if position < cumulative_lengths[i + 1]:
                    offset = position - cumulative_lengths[i]
                    return start + offset
            return -1

        primers = design_primers(sequence, normalized_coords, 20)

        results = []
        for idx, primer_set in enumerate(primers):
            results.append({
                'Pair': idx + 1,
                'Product Size (bp)': primer_set['amplicon_size'],
                'F Pr.': primer_set['left_primer']['sequence'],
                'F Len. (bp)': primer_set['left_primer']['length'],
                'F Pos.': primer_set['left_primer']['position'],
                'F Pos. Abs.': primer_set['left_primer']['position_abs'],
                'F Tm (°C)': primer_set['left_primer']['tm'],
                'F GC%': primer_set['left_primer']['gc_percent'],
                'F Self Compl.': primer_set['left_primer']['self_complementarity'],
                "F Self 3' Compl.": primer_set['left_primer']['self_3prime_complementarity'],
                'R Pr.': primer_set['right_primer']['sequence'],
                'R Len. (bp)': primer_set['right_primer']['length'],
                'R Pos.': primer_set['right_primer']['position'],
                'R Pos. Abs.': primer_set['right_primer']['position_abs'],
                'R Tm (°C)': primer_set['right_primer']['tm'],
                'R GC%': primer_set['right_primer']['gc_percent'],
                'R Self Compl.': primer_set['right_primer']['self_complementarity'],
                "R Self 3' Compl.": primer_set['right_primer']['self_3prime_complementarity'],
                'Product Size Abs (bp)': primer_set['amplicon_size_abs']
            })

        results_df = pd.DataFrame(results)
        print(results_df)

        graphique(normalized_coords, primers)
        graphique(normalized_coords, primers, True)
'''