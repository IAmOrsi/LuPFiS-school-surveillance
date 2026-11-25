import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from matplotlib.patches import Wedge, Circle
from pathlib import Path

# ----------------------------------------
# 1. Configuration
# ----------------------------------------

input_file = "Fig1_data.xlsx"

panels = [
    "Contact",
    "Co-design",
    "Molecular sampling",
    "Lab analysis",
    "Survey",
    "Results workshop"
]

outer_colors = ["white", "#e6e6e6"] * 3

outer_radius = 1
inner_radius = 0.4

node_size = 300
edge_color = "grey"
edge_alpha = 0.8
edge_width = 1.5

# ----------------------------------------
# 2. Function to draw radial background
# ----------------------------------------

def create_pie_background(ax, cloud_color="white"):

    angles = np.linspace(
        np.pi / 6,
        2 * np.pi + np.pi / 6,
        7
    )

    # Outer sectors
    for i, color in enumerate(outer_colors):

        wedge = Wedge(
            (0, 0),
            outer_radius,
            np.degrees(angles[i]),
            np.degrees(angles[i + 1]),
            width=outer_radius - inner_radius,
            facecolor=color,
            edgecolor="none",
            alpha=1,
            zorder=-1
        )

        ax.add_patch(wedge)

    # Central circle
    ax.add_patch(Circle(
        (0, 0),
        inner_radius,
        facecolor="whitesmoke",
        edgecolor="grey",
        zorder=0
    ))

    # Outer circle
    ax.add_patch(Circle(
        (0, 0),
        outer_radius,
        facecolor="none",
        edgecolor="grey",
        zorder=0
    ))

    # White margin
    ax.add_patch(Circle(
        (0, 0),
        outer_radius * 1.1,
        facecolor="white",
        edgecolor="none",
        zorder=-2
    ))

    # External sectors
    for i in [0, 2, 3, 4, 5]:

        cloud = Wedge(
            (0, 0),
            outer_radius * 1.2,
            np.degrees(angles[i]),
            np.degrees(angles[i + 1]),
            width=outer_radius - inner_radius,
            facecolor=cloud_color,
            edgecolor="none",
            alpha=1,
            zorder=-3
        )

        ax.add_patch(cloud)

    ax.set_xlim(-outer_radius * 1.2, outer_radius * 1.2)
    ax.set_ylim(-outer_radius * 1.2, outer_radius * 1.2)

    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


# ----------------------------------------
# 3. Read and plot each panel
# ----------------------------------------

fig, axes = plt.subplots(1, 6, figsize=(23, 5))

with pd.ExcelFile(input_file) as workbook:

    for ax, panel in zip(axes, panels):

        # Read sheet
        data = pd.read_excel(
            workbook,
            sheet_name=panel,
            header=1
        )

        # Extract node table
        nodes = data[
            [
                "Participant",
                "Role",
                "Institute",
                "Group",
                "Radial",
                "X",
                "Y",
                "Active",
                "Alpha",
                "Color"
            ]
        ].dropna(subset=["Participant"]).copy()

        # Extract edge table
        edges = data[
            ["Source", "Target", "Weight"]
        ].dropna(subset=["Source", "Target"]).copy()

        # ----------------------------------------
        # Reconstruct network
        # ----------------------------------------

        G = nx.Graph()

        G.add_nodes_from(nodes["Participant"])

        for _, row in edges.iterrows():

            source = row["Source"]
            target = row["Target"]

            if pd.notna(row["Weight"]):

                G.add_edge(
                    source,
                    target,
                    weight=row["Weight"]
                )

            else:

                G.add_edge(source, target)

        # ----------------------------------------
        # Read visual properties directly
        # from the Excel file
        # ----------------------------------------

        positions = {
            row["Participant"]: (row["X"], row["Y"])
            for _, row in nodes.iterrows()
        }

        node_list = nodes["Participant"].tolist()
        node_colors = nodes["Color"].tolist()
        node_alphas = nodes["Alpha"].tolist()

        # ----------------------------------------
        # Draw background
        # ----------------------------------------

        cloud_color = (
            "#eecc66" if panel == "Survey"
            else "white"
        )

        create_pie_background(
            ax,
            cloud_color=cloud_color
        )

        # ----------------------------------------
        # Draw nodes
        # ----------------------------------------

        nx.draw_networkx_nodes(
            G,
            positions,
            ax=ax,
            nodelist=node_list,
            node_color=node_colors,
            node_size=node_size,
            alpha=node_alphas,
            linewidths=0,
            edgecolors="black"
        )

        # ----------------------------------------
        # Draw edges
        # ----------------------------------------

        nx.draw_networkx_edges(
            G,
            positions,
            ax=ax,
            edge_color=edge_color,
            alpha=edge_alpha,
            width=edge_width
        )

        ax.set_title(panel, fontsize=22)

        # Report network statistics
        n_connections = (
            G.number_of_edges()
            - nx.number_of_selfloops(G)
        )

        print(
            f"{panel}: "
            f"{G.number_of_nodes()} nodes, "
            f"{n_connections} connections"
        )


# ----------------------------------------
# 4. Save figure
# ----------------------------------------

plt.tight_layout()

Path("plots").mkdir(parents=True, exist_ok=True)

plt.savefig("plots/fig1_draft1.pdf")
plt.savefig("plots/fig1_draft1.png", dpi=300)

plt.show()