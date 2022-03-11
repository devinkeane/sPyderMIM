# last rev: 03-03-06


# Import libraries
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import networkx as nx

# Importing fonts --> matplotlib font manager for graph output
import matplotlib.font_manager as font_manager
# Add every font at the specified location
font_dir = ['./fonts/copperplate']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

# Set font family globally
font_manager.rcParams['font.family'] = 'copperplate'


# ------------------------------------------------------------------------------------------------------
# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Apply graph theory to your network table! * ")
parser.add_argument('-m', '--mode', type=str, required=True , help='\"geno\" (find protein interactor overlap) or \"pheno\" (find OMIM phenotypic overlap)')
parser.add_argument('-i', '--input', type=str, required=True ,help='<INPUT_FILENAME.csv>  (Input table)')
parser.add_argument('-l', '--labels', required=False, type=str, default='none', help='Labels --> arguments: \"all\" or \"none\"')
parser.add_argument('-o', '--output', type=str, required=True , help='<OUTPUT_FILENAME.png>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
labels = args.labels
mode = args.mode
# ---------------------------------------------------------------------------
# TITLE SCREEN / LOGO |
# --------------------+

logo = """

O---o    ___|                       _ \   |
O---o   |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
 O-o    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
  O    \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
 o-O   ______________________________________________________________|---------+
o---O   High Performance Computing Genomic Network Analysis    |  Version 4.1  |    ✧ - ･ﾟ*
O---o                               +------------------------------------------+ 
 O-o                     (✿◠‿◠)     |  (c) 2022-01-27 Devin Keane              |
  O                                 |  Feltus Lab                              |◉‿◉)つ
 o-O                                |  Department of Genetics and Biochemistry |
o---O                               |  Clemson University                      | 
O---o                       .'✧     |                                          |
                                    |  Last rev: 2022-03-10                    |
                                    +------------------------------------------+
                         , ⌒ *: ﾟ･✧* ･ﾟ✧ - *                      ─=≡Σ((( つ◕ل͜◕)つ
    ╰( ͡° ͜ʖ ͡° )つ──☆*:・^'
"""
print(logo)

# ---------------------------------------------------------------------------
# GRAPHING THE DATA |
# ------------------+
if mode == 'pheno':
    gpn = pd.read_csv(input)
    print('Processing your input table:')
    print()
    print(gpn.drop(columns='Unnamed: 0'))

    # Create a NetworkX object called "G" where 'Superphenotype' is the source node
    # and 'Node_name' is the target node.
    G = nx.from_pandas_edgelist(gpn,source = 'Superphenotype', target = 'Node_name')

    """
    lst_overlap = gpn[gpn['Node_type'] == 'associated_and_overlapping_genes']['Node_name'].reset_index(drop=True)
    lst_neighbors = gpn[gpn['Node_type'] == 'Intact_first_neighbors']['Node_name'].reset_index(drop=True)
    lst_new = []
    lst_overlap_2 = []
    lst_neighbors_2 = []
    
    
    for i in range(len(lst_overlap)):
        lst_new += [lst_overlap[i]]
        lst_overlap_2 += [lst_overlap[i]]
        lst_new += [lst_neighbors[i]]
        lst_neighbors_2 += [lst_neighbors[i]]
    """

    # Draw a graph with G using a color map that distinguishes between genes and phenotypes
    plt.figure(figsize=(100,100))
    plt.tight_layout()

    """
    # This portion of the code may be deprecated or altered in the future
    # --------------------------------------------------------------------
    for node in G:
        if node in lst_overlap_2:
            color_map.append('yellow')
        elif node in lst_neighbors_2:
            color_map.append('red')
        else:
            color_map.append('green')
    # --------------------------------------------------------------------
    """
    pos = nx.kamada_kawai_layout(G)
    # pos= nx.spring_layout(G)
    nx.draw(G, node_color='green', node_size=300, pos=pos,with_labels=False)



    #color_map = []

    superphenotype_labels = {}
    Node_name_labels = {}
    neighbor_labels = {}
    gene_labels = {}


    if labels == 'all' or labels == 'subtype':
        for idx, node in enumerate(G.nodes()):
            if node in gpn['Superphenotype'].unique():
                superphenotype_labels[node] = node
            if node in gpn['Node_name'].unique():
                Node_name_labels[node] = node

        bbox = dict(fc="blue", ec="black", boxstyle="square", lw=2)

        nx.draw_networkx_labels(G, pos, labels=superphenotype_labels, font_size=14, font_color='white', font_family='copperplate',bbox=bbox)
        nx.draw_networkx_labels(G, pos, labels=Node_name_labels, font_size=14, font_color='black', font_family='copperplate')
    """
    if labels == 'all' or labels == 'overlapping':
        for idx, node in enumerate(G.nodes()):
            if node in lst_overlap_2:
                gene_labels[node] = node

        bbox2 = dict(fc="yellow", ec="black", boxstyle="circle", lw=2)
        nx.draw_networkx_labels(G, pos, labels=gene_labels, font_size=14, font_color='black', font_family='copperplate',bbox=bbox2)

    if labels == 'all' or labels == 'interactors':
        for idx, node in enumerate(G.nodes()):
            if node in lst_neighbors_2:
                neighbor_labels[node] = node

        bbox3 = dict(fc="red", ec="black", boxstyle="sawtooth", lw=2)
        nx.draw_networkx_labels(G, pos, labels=neighbor_labels, font_size=14, font_color='black', font_family='copperplate',bbox=bbox3)

    """
# ---------------------------------------------------------------------------


elif mode == 'geno':
    protein_df = pd.read_csv(input)
    print('Processing your input table:')
    print()
    print(protein_df.drop(columns='Unnamed: 0'))

    # Create a NetworkX object called "G" where 'moleculeA' is the source node
    # and 'moleculeB' is the target node.
    G = nx.from_pandas_edgelist(protein_df, source='moleculeA', target='moleculeB',edge_attr='intactMiscore', create_using=nx.MultiGraph())
    plt.figure(figsize=(100, 100))
    plt.tight_layout()
    pos = nx.kamada_kawai_layout(G)


    if labels == 'all' or labels == 'protein_interactions':
        nx.draw(G, font_color='red', node_color='lightblue', node_size=300, pos=pos, with_labels=True)
    else:
        nx.draw(G, font_color='red', node_color='lightblue', node_size=300, pos=pos, with_labels=False)



# ---------------------------------------------------------------------------
# Save output to files |
# ---------------------+

# Name the graph output file based on the input argument for the file name.
# Append '.png' to the filename and save the figure as that filename.
graph_output_name = output.split('.')[0]
graph_output_name += '.png'
plt.savefig(graph_output_name)

# ---------------------------------------------------------------------------
# Print logo and output message |
# ------------------------------+
num_nodes = G.number_of_nodes()


print()


print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print('     ...Your network graph was saved as \"',output,'\" with ', num_nodes,' total nodes.')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()

spiderweb_ascii_art = """
                                                                 \   ,'|`-.   /
                                                                  \,' _|_  ','
                                                                  /'.' | `,' \\
                                                              -._/_/_'.|,'_\__\_,-
                                                                 | | ,-*." |  |
                                                              ___|,+' /|\`.|  |
             [  G e n o  ]                                       \  \/ | \/`. |___
                                                                  \ /`.|,'\  /
                                                                   Y.  |   \/
                                                                   | `.|_,'
                                                                   |
                                                                   |
                  Network spun.                                 __ |
                                                                __\|,-
                                                                ,-`=--.       
                                                                 /=8\
"""

spiderweb_ascii_art2 = """
                                                                 \   ,'|`-.   /
                                                                  \,' _|_  ','
                                                                  /'.' | `,' \\
                                                              -._/_/_'.|,'_\__\_,-
                                                                 | | ,-*." |  |
                                                              ___|,+' /|\`.|  |
                                    [  P h e n o  ]              \  \/ | \/`. |___
                                                                  \ /`.|,'\  /
                                                                   Y.  |   \/
                                                                   | `.|_,'
                                                                   |
                                                                   |
                  Network spun.                                 __ |
                                                                __\|,-
                                                                ,-`=--.       
                                                                 /=8\
"""

if mode == 'geno':
    print(spiderweb_ascii_art)
if mode == 'pheno':
    print(spiderweb_ascii_art2)
print()