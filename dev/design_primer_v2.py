import primer3
import altair as alt
import pandas as pd


def design_primers(sequence, exons):
    simplified_sequence = "".join(sequence[start:end] for start, end in exons)

    # Paramètres pour Primer3
    primer3_input = {
        'SEQUENCE_ID': 'primer_in_exons',
        'SEQUENCE_TEMPLATE': simplified_sequence,
    }

    primer3_params = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_NUM_RETURN': 2,
        'PRIMER_PRODUCT_SIZE_RANGE': [[80, 250]],  # Adapter à la taille maximale
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_SELF_ANY': 8.0,
        'PRIMER_SELF_END': 3.0,
    }

    primers = []

    exon_lengths = [end - start for start, end in exons]
    cumulative_lengths = [0] + list(cumsum(exon_lengths))

    for i in range(len(exons)):
        for j in range(i + 1, len(exons)):
            exon1_start, exon1_end = exons[i]
            exon2_start, exon2_end = exons[j]

            simplified_start1 = cumulative_lengths[i]
            simplified_end1 = simplified_start1 + (exon1_end - exon1_start)
            simplified_start2 = cumulative_lengths[j]
            simplified_end2 = simplified_start2 + (exon2_end - exon2_start)

            product_size = simplified_end2 - simplified_start1

            if 80 <= product_size <= 250:
                primer3_input['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [
                    simplified_start1, simplified_end1 - simplified_start1,
                    simplified_start2, simplified_end2 - simplified_start2
                ]

                primer_results = primer3.bindings.design_primers(primer3_input, primer3_params)

                if 'PRIMER_PAIR_NUM_RETURNED' in primer_results and primer_results['PRIMER_PAIR_NUM_RETURNED'] > 0:
                    for k in range(primer_results['PRIMER_PAIR_NUM_RETURNED']):
                        left_key = f'PRIMER_LEFT_{k}_SEQUENCE'
                        right_key = f'PRIMER_RIGHT_{k}_SEQUENCE'

                        if left_key in primer_results and right_key in primer_results:
                            left_position = primer_results.get(f'PRIMER_LEFT_{k}')[0]
                            right_position = primer_results.get(f'PRIMER_RIGHT_{k}')[0]

                            left_absolute = convert_to_absolute(left_position, exons, cumulative_lengths)
                            right_absolute = convert_to_absolute(right_position, exons, cumulative_lengths)

                            amplicon_size = right_position - left_position + 1
                            amplicon_size_abs = right_absolute - left_absolute + 1

                            primers.append({
                                'left_primer': {
                                    'sequence': primer_results.get(left_key, 'N/A'),
                                    'position': left_position,
                                    'position_abs': left_absolute,
                                    'tm': primer_results.get(f'PRIMER_LEFT_{k}_TM', 'N/A'),
                                    'gc_percent': primer_results.get(f'PRIMER_LEFT_{k}_GC_PERCENT', 'N/A'),
                                },
                                'right_primer': {
                                    'sequence': primer_results.get(right_key, 'N/A'),
                                    'position': right_position,
                                    'position_abs': right_absolute,
                                    'tm': primer_results.get(f'PRIMER_RIGHT_{k}_TM', 'N/A'),
                                    'gc_percent': primer_results.get(f'PRIMER_RIGHT_{k}_GC_PERCENT', 'N/A'),
                                },
                                'amplicon_size': amplicon_size,
                                'amplicon_size_abs': amplicon_size_abs,

                            })

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


# Exemple d'utilisation
sequence = "gcgcgatcacattccgaactagcaagcggtcgtccagagacggtatcacccatacgtcttgacattacaagttctacctcgagcacacgtattttcacgaaaaacaggattatccaagtctgtcgcactgaacaagtcccccatcagtcctaagtgcgctaccctaccgtcatggggaaccgaacaatcgttggattgggaaagcgagcaaacgtagttccttgtcgagagtgaatttcatttctgatcgcacgagttgcgtaacaggccagatggtcttttatgaaggcggacggccagggttgcatggctctccgtgaatccggagaagccttctaataattatcctcttgactaatgtcggccacaaagtggttgtggttggtatctttttcgtcagttagccatgccagtcggaacgagcagcagggtattctaccgaaaaccggcccagctcaagttgcctaatctgaaccctatatcgggtttagtcgaaatgtagcatgaatcgtggactggcgctgcattccgggacgcgttccatcgcgtgcagccgaaacgctgagtgagaaagtttctgtccgctttcctggctctccgcgacccttctgtctgctaactcctgccgcgacccgacaacggaactgccttcaggtgccccccccagaggaagggtcttggac"

# Liste des exons
exons = [(0, 150), (250, 350), (425, 750)]
# Conception des amorces
primers = design_primers(sequence, exons)

# Affichage des résultats
for idx, primer_set in enumerate(primers):
    print(f"\nJeu de primers {idx + 1}:")
    print(f"  Forward Primer: {primer_set['left_primer']}")
    print(f"  Reverse Primer: {primer_set['right_primer']}")
    print(f"  Amplicon Size: {primer_set['amplicon_size']} bp")
    print(f"  Amplicon Size Abs: {primer_set['amplicon_size_abs']} bp")

exons_data = pd.DataFrame(exons, columns=['start', 'end'])
primers_data = [
    {'type': 'Left Primer', 'position_abs': p['left_primer']['position_abs'], 'label': f'Left {i}'}
    for i, p in enumerate(primers)
] + [
    {'type': 'Right Primer', 'position_abs': p['right_primer']['position_abs'], 'label': f'Right {i}'}
    for i, p in enumerate(primers)
]
primers_data = pd.DataFrame(primers_data)

# Graphiques Altair
exons_chart = alt.Chart(exons_data).mark_rect(stroke='darkblue', strokeWidth=2).encode(
    x=alt.X('start:Q', title='Position dans la Séquence'),
    x2='end:Q',
    color=alt.value('lightblue'),
    tooltip=['start', 'end']
).interactive()

primers_chart = alt.Chart(primers_data).mark_point(filled=True, size=100).encode(
    x=alt.X('position_abs:Q'),
    y=alt.value(10),
    color=alt.Color('type:N', scale=alt.Scale(domain=['Left Primer', 'Right Primer'], range=['red', 'green'])),
    tooltip=['label', 'position_abs']
)

annotations = alt.Chart(primers_data).mark_text(
    align='left', dx=10, dy=-10, fontSize=10, color='black'
).encode(
    x='position_abs:Q',
    y=alt.value(0),
    text='label:N'
)

final_chart = exons_chart + primers_chart + annotations

# Sauvegarder ou afficher
final_chart.save('sequence_visualization.html')
print("Graphique sauvegardé sous 'sequence_visualization.html'. Ouvrez ce fichier pour le visualiser.")

# Préparation des données pour les exons
adjusted_exons = []
intron_cumulative_length = 0
for i, (start, end) in enumerate(exons):
    # Si ce n'est pas le premier exon, calculer l'intron précédent
    if i > 0:
        intron_cumulative_length += exons[i][0] - exons[i - 1][1]  # Distance entre exon[i] et exon[i-1]

    # Ajuster les coordonnées de l'exon en fonction de l'intron cumulé
    adjusted_start = start - intron_cumulative_length
    adjusted_end = end - intron_cumulative_length
    adjusted_exons.append((adjusted_start, adjusted_end))

exons_data_adj = pd.DataFrame(adjusted_exons, columns=['start', 'end'])

# Création d'un champ 'pair' pour les primers, en assignant un même jeu à chaque couple de primers (gauche + droite)
primers_data_adj = [
    {'pair': f'Pair {i + 1}', 'type': 'Left Primer', 'position': p['left_primer']['position'], 'label': f'Left {i + 1}', 'y': (i + 1) * 10}
    for i, p in enumerate(primers)
] + [
    {'pair': f'Pair {i + 1}', 'type': 'Right Primer', 'position': p['right_primer']['position'], 'label': f'Right {i + 1}', 'y': (i + 1) * 10}
    for i, p in enumerate(primers)
]

primers_data_adj = pd.DataFrame(primers_data_adj)

# Graphiques Altair

# Exons chart avec y dynamique
exons_chart_adj = alt.Chart(exons_data_adj).mark_rect(stroke='darkblue', strokeWidth=2).encode(
    x=alt.X('start:Q', title='Position dans la Séquence'),
    x2='end:Q',
    color=alt.value('lightblue'),
    tooltip=['start', 'end']
).interactive()

# Primers chart, maintenant colorié par jeu de primers
primers_chart_adj = alt.Chart(primers_data_adj).mark_point(filled=True, size=100).encode(
    x=alt.X('position:Q'),
    y=alt.Y('y:Q', title=None),  # Utilise la colonne 'y' pour espacer les points
    color=alt.Color('pair:N', scale=alt.Scale(domain=[f'Pair {i + 1}' for i in range(len(primers))]), title='Primers pairs'),  # Couleur par jeu
    tooltip=['label', 'position']
)

# Annotations des primers
annotations_adj = alt.Chart(primers_data_adj).mark_text(
    align='left', dx=10, dy=-10, fontSize=10, color='black'
).encode(
    x='position:Q',
    y='y:Q',
    text='label:N'
)

# Combiner Exons et Primers
final_chart_adj = exons_chart_adj + primers_chart_adj + annotations_adj

# Sauvegarder ou afficher
final_chart_adj.save('sequence_visualization_adjusted.html')
print("Graphique sauvegardé sous 'sequence_visualization_adjusted.html'. Ouvrez ce fichier pour le visualiser.")
