from typing import List, Optional, Union

import altair as alt
import pandas as pd

from .components import create_horizontal_bar, create_matrix_view, create_vertical_bar
from .config import upsetaltair_top_level_configuration
from .preprocessing import preprocess_data
from .transforms import create_base_chart


class UpSetChart:
    """A wrapper class for UpSet plots."""

    def __init__(self, chart, data, sets):
        """Initialize the UpSetChart.

        Parameters
        ----------
        chart : alt.Chart
            The base Altair chart
        data : pd.DataFrame
            The input data
        sets : list
            List of set names
        """
        self.chart = chart
        self.data = data
        self.sets = sets

    def save(self, filename, format="png"):
        if format not in['png', 'svg', 'pdf', 'html', 'json', 'vega']:
            raise ValueError("Format not supported, must specify file format: ['png', 'svg', 'pdf', 'html', 'json', 'vega']")
        """Save the chart to a file."""
        self.chart.save(filename, format=format)

    def properties(self, **kwargs):
        """Update chart properties."""
        self.chart = self.chart.properties(**kwargs)
        return self

    def configure_axis(self, **kwargs):
        """Configure chart axes."""
        self.chart = self.chart.configure_axis(**kwargs)
        return self

    def configure_legend(self, **kwargs):
        """Configure chart legend."""
        self.chart = self.chart.configure_legend(**kwargs)
        return self

    def to_dict(self):
        """Convert the chart to a dictionary representation.

        Returns
        -------
        dict
            The Vega-Lite specification as a Python dictionary
        """
        return self.chart.to_dict()

    def __getattr__(self, name):
        """Delegate unknown attributes to the underlying chart."""
        return getattr(self.chart, name)


def UpSetAltair(
    data: pd.DataFrame,
    sets: List[str],
    *,
    title: str = "",
    subtitle: Union[str, List[str]] = "",
    abbre: Optional[List[str]] = None,
    sort_by: str = "frequency",
    sort_order: str = "ascending",
    width: int = 1200,
    height: int = 700,
    height_ratio: float = 0.6,
    horizontal_bar_chart_width: Optional[int] = None,
    color_range: List[str] = [
        "#55A8DB",
        "#3070B5",
        "#30363F",
        "#F1AD60",
        "#DF6234",
        "#BDC6CA",
    ],
    highlight_color: str = "#EA4667",
    glyph_size: int = 100,  # Reduced from 200
    set_label_bg_size: int = 500,  # Reduced from 1000
    line_connection_size: int = 1,  # Reduced from 2
    horizontal_bar_size: int = 20,
    vertical_bar_label_size: int = 16,
    # vertical_bar_padding: int = 20,
    theme: Optional[str] = None,
) -> UpSetChart:
    """Generate interactive UpSet plots using Altair. [Lex et al., 2014]_

    UpSet plots are used to visualize set intersections in a more scalable way than Venn diagrams.
    This implementation provides interactive features like hover highlighting and legend filtering.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data where each column represents a set and contains binary values (0 or 1).
        Each row represents an element, and the columns indicate set membership.
    sets : list of str
        Names of the sets to visualize (must correspond to column names in data).
    title : str, default ""
        Title of the plot.
    subtitle : str or list of str, default ""
        Subtitle(s) of the plot. Can be a single string or list of strings for multiple lines.
    abbre : list of str, optional
        Abbreviations for set names (must have same length as sets).
    sort_by : {"frequency", "degree"}, default "frequency"
        Method to sort the intersections:
        - "frequency": sort by intersection size
        - "degree": sort by number of sets in intersection
    sort_order : {"ascending", "descending"}, default "ascending"
        Order of sorting for intersections.
    width : int, default 1200
        Total width of the plot in pixels.
    height : int, default 700
        Total height of the plot in pixels.
    height_ratio : float, default 0.6
        Ratio of vertical bar chart height to total height (between 0 and 1).
    horizontal_bar_chart_width : int, default 300
        Width of the horizontal bar chart in pixels.
    color_range : list of str
        List of colors for the sets. Defaults to a colorblind-friendly palette.
    highlight_color : str, default "#EA4667"
        Color used for highlighting on hover.
    glyph_size : int, default 200
        Size of the matrix glyphs in pixels.
    set_label_bg_size : int, default 1000
        Size of the set label background circles.
    line_connection_size : int, default 2
        Thickness of connecting lines in pixels.
    horizontal_bar_size : int, default 20
        Height of horizontal bars in pixels.
    vertical_bar_label_size : int, default 16
        Font size of vertical bar labels.
    vertical_bar_padding : int, default 20
        Padding between vertical bars.
    theme : str, optional
        Altair theme to use. If None, uses the current default theme.

    Returns
    -------
    altair.Chart
        An Altair chart object representing the UpSet plot.

    Examples
    --------
    >>> import altair_upset as au
    >>> import pandas as pd
    >>> data = pd.DataFrame({
    ...     'set1': [1, 0, 1],
    ...     'set2': [1, 1, 0],
    ...     'set3': [0, 1, 1]
    ... })
    >>> chart = au.UpSetAltair(
    ...     data=data,
    ...     sets=["set1", "set2", "set3"],
    ...     title="Sample UpSet Plot"
    ... )

    Notes
    -----
    The plot consists of three main components:
    1. A matrix view showing set intersections
    2. A vertical bar chart showing intersection sizes
    3. A horizontal bar chart showing set sizes

    References
    ----------
    .. [Lex et al., 2014] Alexander Lex, Nils Gehlenborg, Hendrik Strobelt, Romain Vuillemot, Hanspeter Pfister.
                UpSet: Visualization of Intersecting Sets
                IEEE transactions on visualization and computer graphics, 20(12), 1983-1992.
    """
    # Input validation
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    if not isinstance(sets, list) or not all(isinstance(s, str) for s in sets):
        raise TypeError("sets must be a list of strings")
    if not all(s in data.columns for s in sets):
        raise ValueError("all sets must be columns in data")
    if not all(data[s].isin([0, 1]).all() for s in sets):
        raise ValueError("all set columns must contain only 0s and 1s")
    if height_ratio <= 0 or height_ratio >= 1:
        raise ValueError("height_ratio must be between 0 and 1")
    if sort_by not in ["frequency", "degree"]:
        raise ValueError("sort_by must be either 'frequency' or 'degree'")
    if sort_order not in ["ascending", "descending"]:
        raise ValueError("sort_order must be either 'ascending' or 'descending'")
    if abbre is not None and len(sets) != len(abbre):
        raise ValueError("if provided, abbre must have the same length as sets")

    # Apply theme if specified
    if theme is not None:
        alt.themes.enable(theme)

    # Preprocess data
    data, set_to_abbre, set_to_order, abbre = preprocess_data(
        data, sets, abbre, sort_order
    )

    # Setup selections for interactivity
    legend_selection = alt.selection_point(fields=["set"], bind="legend")
    color_selection = alt.selection_point(fields=["intersection_id"], on="mouseover")
    opacity_selection = alt.selection_point(fields=["intersection_id"])

    # Calculate dimensions
    if horizontal_bar_chart_width is None:
        horizontal_bar_chart_width = int(width * 0.15)  # Make it 25% of total width
    vertical_bar_chart_height = height * height_ratio
    matrix_height = (height - vertical_bar_chart_height) * 0.8  # Reduce height to tighten spacing
    matrix_width = width - horizontal_bar_chart_width

    # Future update
    num_intersections = max(1, len(data["intersection_id"].unique().tolist()))
    # vertical_bar_size = min(30, max(0, width / num_intersections - vertical_bar_padding))

    # Automatic padding
    vertical_bar_size = min(30, (matrix_width / num_intersections) - 5)

    # vertical_bar_size = min(
    #     30,
    #     width / len(data["intersection_id"].unique().tolist()) - vertical_bar_padding,
    # )

    # Setup styles
    main_color = "#3A3A3A"
    brush_color = alt.condition(
        ~color_selection, alt.value(main_color), alt.value(highlight_color)
    )
    is_show_horizontal_bar_label_bg = len(abbre[0]) <= 2 if abbre else True
    horizontal_bar_label_bg_color = (
        "white" if is_show_horizontal_bar_label_bg else "black"
    )
    x_sort = alt.Sort(
        field="count" if sort_by == "frequency" else "degree", order=sort_order
    )
    tooltip = [
        alt.Tooltip("count:Q", title="Cardinality"),
        alt.Tooltip("degree:Q", title="Degree"),
        # alt.Tooltip("sets_graph:N", title="Groups"),
    ]

    # Create base chart
    base = create_base_chart(data, sets, legend_selection, set_to_abbre, set_to_order)

    # Create components
    vertical_bar, vertical_bar_text = create_vertical_bar(
        base,
        matrix_width,
        vertical_bar_chart_height,
        main_color,
        vertical_bar_size,
        brush_color,
        x_sort,
        tooltip,
        vertical_bar_label_size,
    )
    vertical_bar_chart = (
        (vertical_bar + vertical_bar_text)
        .add_params(color_selection)
        .properties(width=matrix_width, height=vertical_bar_chart_height)
    )

    circle_bg, rect_bg, circle, line_connection = create_matrix_view(
        vertical_bar,
        matrix_height,
        glyph_size,
        x_sort,
        brush_color,
        line_connection_size,
        main_color,
    )
    matrix_view = (
        (circle + rect_bg + circle_bg + line_connection + circle)
        .add_params(color_selection)
        .properties(width=matrix_width)
    )

    horizontal_bar_label_bg, horizontal_bar_label, horizontal_bar = (
        create_horizontal_bar(
            base,
            set_label_bg_size,
            sets,
            color_range,
            is_show_horizontal_bar_label_bg,
            horizontal_bar_label_bg_color,
            horizontal_bar_size,
            horizontal_bar_chart_width,
        )
    )
    horizontal_bar_axis = (
        (horizontal_bar_label_bg + horizontal_bar_label)
        if is_show_horizontal_bar_label_bg
        else horizontal_bar_label
    ).properties(width=horizontal_bar_chart_width)

    # Combine components
    upsetaltair = alt.vconcat(
        vertical_bar_chart,
        alt.hconcat(
            matrix_view,
            horizontal_bar_axis,
            horizontal_bar.properties(width=horizontal_bar_chart_width),
            spacing=0,  # Minimize spacing between components
        ).resolve_scale(x="shared", y="shared"),  # X shared also
        spacing=5,
    ).add_params(legend_selection)

    # Apply configuration
    chart = upsetaltair_top_level_configuration(
        upsetaltair, legend_orient="top", legend_symbol_size=set_label_bg_size / 2.0
    ).properties(
        title={
            "text": title,
            "subtitle": subtitle,
            "fontSize": 20,
            "fontWeight": 500,
            "subtitleColor": main_color,
            "subtitleFontSize": 14,
        }
    )

    return UpSetChart(chart, data, sets)
