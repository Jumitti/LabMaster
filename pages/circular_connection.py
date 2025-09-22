import re
from io import BytesIO

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from utils.page_config import page_config

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


def draw_text(ax, x, y, text, width, height, max_fontsize=10, min_fontsize=6, zorder=4):

    renderer = ax.figure.canvas.get_renderer()
    fontsize = max_fontsize

    words = re.split(r'\s+|-', text)

    while fontsize >= min_fontsize:
        lines = []
        current_line = ""
        for w in words:
            trial_line = f"{current_line} {w}".strip() if current_line else w

            text_obj_temp = ax.text(0, 0, trial_line, fontsize=fontsize,
                                    fontweight='bold', zorder=zorder)
            bbox = text_obj_temp.get_window_extent(renderer=renderer)
            text_obj_temp.remove()

            if bbox.width <= width:
                current_line = trial_line
            else:
                if current_line:
                    lines.append(current_line)
                current_line = w
        if current_line:
            lines.append(current_line)

        wrapped_text = "\n".join(lines)

        text_obj_temp = ax.text(0, 0, wrapped_text, fontsize=fontsize,
                                fontweight='bold', zorder=zorder, linespacing=1.1)
        bbox = text_obj_temp.get_window_extent(renderer=renderer)
        text_obj_temp.remove()
        if bbox.height <= height:
            break
        fontsize -= 1

    ax.text(x, y, wrapped_text, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', zorder=zorder, linespacing=1.1)


page_config()
st.title("⭕ Circular graph (Chord Style)")

if "elements" not in st.session_state:
    st.session_state["elements"] = []
if "edges" not in st.session_state:
    st.session_state["edges"] = {}

new_element = st.text_input("➕ Add elements :", "")
if st.button("Add element") and new_element and new_element not in st.session_state["elements"]:
    st.session_state["elements"].append(new_element)
    st.session_state["edges"][new_element] = []

if st.session_state["elements"]:
    to_delete = st.multiselect("❌ Delete elements :", options=st.session_state["elements"])
    if st.button("Delete elements") and to_delete:
        for el in to_delete:
            st.session_state["elements"].remove(el)
            st.session_state["edges"].pop(el, None)
        for el, conns in st.session_state["edges"].items():
            st.session_state["edges"][el] = [c for c in conns if c not in to_delete]


if st.session_state["elements"]:
    st.markdown("### 🎯 Elements :")
    st.write(", ".join(st.session_state["elements"]))

for el in st.session_state["elements"]:
    others = [o for o in st.session_state["elements"] if o != el]
    current_connections = st.session_state["edges"].get(el, [])

    selected = st.multiselect(
        f"🔗 Connected with {el} :",
        options=others,
        default=current_connections,
        key=f"ms_{el}"
    )

    st.session_state["edges"][el] = selected
    for s in selected:
        if el not in st.session_state["edges"][s]:
            st.session_state["edges"][s].append(el)
    for s in list(st.session_state["edges"].keys()):
        if s != el and el in st.session_state["edges"][s] and s not in selected:
            st.session_state["edges"][s].remove(el)

if st.session_state["elements"]:
    cmap_choice = st.selectbox("🎨 Colors :", list(cmap_options.keys()), index=0)
    cmap = plt.get_cmap(cmap_options[cmap_choice])
    colors = [cmap(i / max(1, len(st.session_state["elements"]) - 1)) for i in range(len(st.session_state["elements"]))]

    radius = 5.0
    circle_color = st.color_picker("🎯 Circle color", "#000000")
    circle_lw = st.slider("Circle width", min_value=0.5, max_value=6.0, value=2.0, step=0.5)
    curvature = st.slider("📐 Innerness", 0.0, 1.0, 0.4, 0.05)

    elements = st.session_state["elements"]
    edges = []
    for src, dsts in st.session_state["edges"].items():
        for dst in dsts:
            if (dst, src) not in edges:
                edges.append((src, dst))

    n = len(elements)
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    positions = {el: (radius*np.cos(t), radius*np.sin(t)) for el, t in zip(elements, theta)}

    fig, ax = plt.subplots(figsize=(7, 7), dpi=300)

    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    ax.set_aspect("equal")
    ax.axis("off")

    circle_patch = plt.Circle((0, 0), radius, edgecolor=circle_color, facecolor="none", lw=circle_lw, zorder=0)
    ax.add_patch(circle_patch)

    margin = radius * 0.2
    ax.set_xlim(-radius - margin, radius + margin)
    ax.set_ylim(-radius - margin, radius + margin)

    shape_option = st.radio("🔹 Elements form :", ["Circle", "Rectangle"])
    node_size = st.slider("Elements size", min_value=0.5, max_value=2.5, value=1.2, step=0.1)

    for (el, (x, y)), c in zip(positions.items(), colors):
        width = height = 1.0 * node_size

        if shape_option == "Circle":
            ax.scatter(x, y, s=1200 * node_size, color=c, edgecolors="black", zorder=3)
        else:
            rect = mpatches.FancyBboxPatch(
                (x - width / 2, y - height / 2),
                width, height,
                boxstyle="round,pad=0.1",
                facecolor=c,
                edgecolor="black",
                zorder=3
            )
            ax.add_patch(rect)

        draw_text(ax, x, y, el, width, height, max_fontsize=14, min_fontsize=6)

    for src, dst in edges:
        if src in positions and dst in positions:
            x1, y1 = positions[src]
            x2, y2 = positions[dst]
            xm, ym = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            control = (xm * curvature, ym * curvature)
            path = mpatches.Path(
                [(x1, y1), control, (x2, y2)],
                [mpatches.Path.MOVETO, mpatches.Path.CURVE3, mpatches.Path.CURVE3]
            )
            patch = mpatches.PathPatch(path, facecolor="none", edgecolor="gray", lw=1.4, alpha=0.8, zorder=2)
            ax.add_patch(patch)

    st.pyplot(fig)

    buffer = BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight", dpi=300, transparent=True)
    buffer.seek(0)

    st.download_button(
        label="💾 Download (.png)",
        data=buffer,
        file_name="circular_graph_tags.png",
        mime="image/png"
    )
