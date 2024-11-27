import requests
import xml.etree.ElementTree as ET
import time
import random
from Bio.Seq import Seq
from Bio import Entrez
import primer3
from dev.tfinder import NCBIdna
import altair as alt
import pandas as pd
from tqdm import tqdm


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
                # Tâche principale
                'PRIMER_TASK': 'generic',  # Génération d'amorces génériques

                # Paramètres pour les amorces (taille et contenu)
                'PRIMER_OPT_SIZE': 20,  # Taille optimale de l'amorce
                'PRIMER_MIN_SIZE': 18,  # Taille minimale de l'amorce
                'PRIMER_MAX_SIZE': 24,  # Taille maximale de l'amorce
                'PRIMER_OPT_TM': 60.0,  # Température de fusion optimale (°C)
                'PRIMER_MIN_TM': 57.0,  # Température de fusion minimale (°C)
                'PRIMER_MAX_TM': 63.0,  # Température de fusion maximale (°C)
                'PRIMER_MIN_GC': 45.0,  # Pourcentage GC minimum (%)
                'PRIMER_MAX_GC': 55.0,  # Pourcentage GC maximum (%)
                'PRIMER_GC_CLAMP': 0,  # Clampage GC à l'extrémité 3' (nombre minimum de G/C)
                'PRIMER_MAX_POLY_X': 5,  # Nombre maximal de bases répétées (ex : AAAAA)

                # Paramètres pour la stabilité à l'extrémité 3'
                'PRIMER_MAX_END_STABILITY': 9.0,  # Stabilité maximale de l'extrémité 3'

                # Alignement secondaire (Thermodynamic model)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 70.0,  # Mauvais appariement au modèle (paires d'amorces)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH_TMPL': 40.0,  # Mauvais appariement pour une amorce unique
                'PRIMER_MAX_SELF_ANY_TH': 45.0,  # Appariement interne (tous les sites, thermodynamique)
                'PRIMER_MAX_SELF_END_TH': 35.0,  # Appariement interne (extrémité 3', thermodynamique)
                'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,  # Appariement entre amorces (tous les sites, thermodynamique)
                'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,  # Appariement entre amorces (extrémité 3', thermodynamique)
                'PRIMER_MAX_HAIRPIN_TH': 24.0,  # Énergie libre maximale pour les épingles à cheveux

                # Alignement secondaire (Ancien modèle)
                'PRIMER_MAX_TEMPLATE_MISPRIMING': 24.0,  # Mauvais appariement au modèle (paires d'amorces, classique)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TMPL': 12.0,  # Mauvais appariement pour une amorce unique (classique)
                'PRIMER_MAX_SELF_ANY': 8.0,  # Appariement interne (tous les sites, classique)
                'PRIMER_MAX_SELF_END': 3.0,  # Appariement interne (extrémité 3', classique)
                'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,  # Appariement entre amorces (tous les sites, classique)
                'PRIMER_PAIR_MAX_COMPL_END': 3.0,  # Appariement entre amorces (extrémité 3', classique)

                # Recherche des alignements secondaires
                'PRIMER_THERMODYNAMIC_ALIGNMENT': 1,  # Utiliser le modèle thermodynamique
                'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1,
                # Aligner également avec le modèle thermodynamique (peut être lent)

                # Paramètres généraux pour les paires
                'PRIMER_NUM_RETURN': 0,  # Nombre maximal de paires retournées
                'PRIMER_PRODUCT_SIZE_RANGE': [[80, 250]],  # Plage de tailles des produits

                # Activer la sélection des oligos d’hybridation internes
                'PRIMER_PICK_INTERNAL_OLIGO': 0,  # 1 pour activer, 0 pour désactiver

                # Paramètres de taille pour les oligos internes
                'PRIMER_INTERNAL_MIN_SIZE': 18,  # Taille minimale
                'PRIMER_INTERNAL_OPT_SIZE': 20,  # Taille optimale
                'PRIMER_INTERNAL_MAX_SIZE': 27,  # Taille maximale

                # Température de fusion (Tm) pour les oligos internes
                'PRIMER_INTERNAL_MIN_TM': 57.0,  # Tm minimale (°C)
                'PRIMER_INTERNAL_OPT_TM': 60.0,  # Tm optimale (°C)
                'PRIMER_INTERNAL_MAX_TM': 63.0,  # Tm maximale (°C)

                # Pourcentage GC pour les oligos internes
                'PRIMER_INTERNAL_MIN_GC': 20.0,  # Pourcentage GC minimal
                'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,  # Pourcentage GC optimal
                'PRIMER_INTERNAL_MAX_GC': 80.0,  # Pourcentage GC maximal

                # Concentration des cations monovalents (ex: Na+)
                'PRIMER_MONOVALENT_CATION_CONC': 50.0,  # En mM (par défaut c'est 50 mM)

                # Concentration des cations divalents (ex: Mg2+)
                'PRIMER_DIVALENT_CATION_CONC': 1.5,  # En mM (par défaut c'est 1.5 mM)

                # Concentration des dNTPs
                'PRIMER_DNTP_CONC': 0.6,  # En mM (par défaut c'est 0.6 mM)

                # Formule de correction du sel (utilisation de la formule de Santa Lucia 1998)
                'PRIMER_SALT_CORRECTION': 1,  # 1 pour utiliser la correction de Santa Lucia 1998

                # Paramètres thermodynamiques (tableau de paramètres pour calculer les Tm)
                'PRIMER_THERMODYNAMIC_PARAMETERS': 'SantaLucia1998',  # Utiliser la table de Santa Lucia 1998

                # Concentration de l'oligonucléotide pour l'initiation de l'hybridation
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