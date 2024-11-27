import altair as alt
import pandas as pd

# Exemple de séquence
sequence = "gcgcgatcacattccgaactagcaagcggtcgtccagagacggtatcacccatacgtcttgacattacaagttctacctcgagcacacgtattttcacgaaaaacaggattatccaagtctgtcgcactgaacaagtcccccatcagtcctaagtgcgctaccctaccgtcatggggaaccgaacaatcgttggattgggaaagcgagcaaacgtagttccttgtcgagagtgaatttcatttctgatcgcacgagttgcgtaacaggccagatggtcttttatgaaggcggacggccagggttgcatggctctccgtgaatccggagaagccttctaataattatcctcttgactaatgtcggccacaaagtggttgtggttggtatctttttcgtcagttagccatgccagtcggaacgagcagcagggtattctaccgaaaaccggcccagctcaagttgcctaatctgaaccctatatcgggtttagtcgaaatgtagcatgaatcgtggactggcgctgcattccgggacgcgttccatcgcgtgcagccgaaacgctgagtgagaaagtttctgtccgctttcctggctctccgcgacccttctgtctgctaactcctgccgcgacccgacaacggaactgccttcaggtgccccccccagaggaagggtcttggac"

# Liste des exons
exons = [(0, 150), (250, 350), (425, 750)]

# Liste simulée de primers générés
primers = [
    {'left_primer': {'position_abs': 0}, 'right_primer': {'position_abs': 150}},
    {'left_primer': {'position_abs': 250}, 'right_primer': {'position_abs': 350}},
    {'left_primer': {'position_abs': 425}, 'right_primer': {'position_abs': 750}},
]

# Préparer les données pour Altair
exons_data = pd.DataFrame(exons, columns=['start', 'end'])
primers_data = [
    {'type': 'Left Primer', 'position': p['left_primer']['position_abs'], 'label': f'Left {i}'}
    for i, p in enumerate(primers)
] + [
    {'type': 'Right Primer', 'position': p['right_primer']['position_abs'], 'label': f'Right {i}'}
    for i, p in enumerate(primers)
]
primers_data = pd.DataFrame(primers_data)

# Graphiques Altair
exons_chart = alt.Chart(exons_data).mark_bar(height=20).encode(
    x=alt.X('start:Q', title='Position dans la Séquence'),
    x2='end:Q',
    color=alt.value('lightblue'),
    tooltip=['start', 'end']
)

primers_chart = alt.Chart(primers_data).mark_point(filled=True, size=100).encode(
    x=alt.X('position_abs:Q'),
    y=alt.value(0),
    color=alt.Color('type:N', scale=alt.Scale(domain=['Left Primer', 'Right Primer'], range=['red', 'green'])),
    tooltip=['label', 'position']
)

annotations = alt.Chart(primers_data).mark_text(
    align='left', dx=10, dy=-10, fontSize=10, color='black'
).encode(
    x='position:Q',
    y=alt.value(0),
    text='label:N'
)

final_chart = exons_chart + primers_chart + annotations

# Sauvegarder ou afficher
final_chart.save('sequence_visualization.html')
print("Graphique sauvegardé sous 'sequence_visualization.html'. Ouvrez ce fichier pour le visualiser.")
