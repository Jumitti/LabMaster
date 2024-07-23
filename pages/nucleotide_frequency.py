import streamlit as st
import pandas as pd
import altair as alt

# Fonction pour calculer les fréquences des nucléotides tous les 5 nucléotides
def calculate_frequencies(sequence, window_size=50):
    frequencies = []
    for i in range(0, len(sequence), window_size):
        window = sequence[i:i + window_size]
        base_counts = {base: window.count(base) for base in "CGAT"}
        total_bases = len(window)
        for base, count in base_counts.items():
            frequencies.append({
                'position': i,
                'base': base,
                'frequency': count / total_bases if total_bases > 0 else 0
            })
    return pd.DataFrame(frequencies)

# Configuration de Streamlit
st.title("Visualisation de Séquence en Ruban")
sequence = st.text_input("Entrez la séquence", "ATGCGCATGCGC")

# Calcul des fréquences
frequencies = calculate_frequencies(sequence)

# Définir les couleurs des bases
base_color = {
    "C": "blue",
    "G": "green",
    "A": "red",
    "T": "orange"
}

# Création du graphique avec Altair
chart = alt.Chart(frequencies).mark_area().encode(
    x=alt.X('position:Q', axis=alt.Axis(title='Position')),
    y=alt.Y('frequency:Q', stack='zero', axis=alt.Axis(title='Fréquence des nucléotides')),
    color=alt.Color('base:N', scale=alt.Scale(domain=list(base_color.keys()), range=list(base_color.values()))),
    order=alt.Order('base:N')
).properties(
    width=800,
    height=400,
    title="Fréquence des Nucléotides tous les 5 Nucléotides"
)

# Affichage du graphique dans Streamlit
st.altair_chart(chart, use_container_width=True)
