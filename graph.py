#                               ,--.
#   ,---. ,--.--. ,--,--. ,---. |  ,---.      ,---.,--. ,--.
#  | .-. ||  .--'' ,-.  || .-. ||  .-.  |    | .-. |\  '  /
#  ' '-' '|  |   \ '-'  || '-' '|  | |  |.--.| '-' ' \   '
#  .`-  / `--'    `--`--'|  |-' `--' `--''--'|  |-'.-'  /
#  `---'                 `--'                `--'  `---'
#                                [  G e n o P h e n o  ]
#
# last rev: 03-03-13
# ------------------------------------------------------------------------------------------------------

# Import libraries
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import networkx as nx
from operator import itemgetter, attrgetter
from datetime import datetime
from networkx.algorithms.flow import shortest_augmenting_path
import itertools
import threading
import time
import sys

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
# ------------------------------------------------------------------------------------------------------
# Importing fonts --> matplotlib font manager for graph output
import matplotlib.font_manager as font_manager
# Add every font at the specified location
font_dir = ['./fonts/copperplate']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

# Set font family globally
font_manager.rcParams['font.family'] = 'copperplate'
# ------------------------------------------------------------------------------------------------------
done = False
# Loading bar animation function
def animate():


    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            break
        sys.stdout.write('\r( ͡°_ʖ ͡°) Calculating... ' + c)
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write('\r( ͡° ͜ʖ ͡°)_/¯ Calculation complete! ✔')


calculation_wait_animation = threading.Thread(target=animate)

# ------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# TITLE SCREEN / LOGO |
# --------------------+

logo = """

O---o    ___|                       _ \   |
O---o   |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
 O-o    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
  O    \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
 o-O   ______________________________________________________________|---------+
o---O   High Performance Computing Genomic Network Analysis    |  Version 5.5  |    ✧ - ･ﾟ*
O---o                               +------------------------------------------+ 
 O-o                     (✿◠‿◠)     |  (c) 2022-01-27 Devin Keane              |
  O                                 |  Feltus Lab                              |◉‿◉)つ
 o-O                                |  Department of Genetics and Biochemistry |
o---O                               |  Clemson University                      | 
O---o                       .'✧     |                                          |
                                    |  Last rev: 2022-04-25                    |
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

    print()
    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'Superphenotype' is the source node
    # and 'Node_name' is the target node.
    G = nx.from_pandas_edgelist(gpn,source = 'Superphenotype', target = 'Node_name')

    sys.stdout.flush()
    sys.stdout.write('\rDefining source and target nodes... ✔')
    print()

    sys.stdout.write('Creating layout...')
    # Draw a graph with G using a color map that distinguishes between genes and phenotypes
    plt.figure(figsize=(100,100))
    sys.stdout.flush()
    plt.tight_layout()





    """
    # This portion of the code may be deprecated or altered in the future
    # --------------------------------------------------------------------
    
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
    sys.stdout.flush()
    sys.stdout.write('\rCreating layout... ✔')
    print()

    sys.stdout.write('Constructing network visualization...')

    nx.draw(G, node_color='green', node_size=300, pos=pos,with_labels=False)


    #color_map = []
    sys.stdout.flush()
    sys.stdout.write('\rConstructing network visualization... ✔')
    print()
    sys.stdout.flush()
    superphenotype_labels = {}
    Node_name_labels = {}
    neighbor_labels = {}
    gene_labels = {}



    if labels == 'all' or labels == 'subtype':
        sys.stdout.write('Creating labels...')
        for idx, node in enumerate(G.nodes()):
            if node in gpn['Superphenotype'].unique():
                superphenotype_labels[node] = node
            if node in gpn['Node_name'].unique():
                Node_name_labels[node] = node

        bbox = dict(fc="blue", ec="black", boxstyle="square", lw=2)



        nx.draw_networkx_labels(G, pos, labels=superphenotype_labels, font_size=14, font_color='white', font_family='copperplate',bbox=bbox)
        nx.draw_networkx_labels(G, pos, labels=Node_name_labels, font_size=14, font_color='black', font_family='copperplate')

        sys.stdout.write('\rCreating labels... ✔')
        sys.stdout.flush()
        print()
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

    print()
    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'moleculeA' is the source node
    # and 'moleculeB' is the target node.
    G = nx.from_pandas_edgelist(protein_df, source='moleculeA', target='moleculeB',edge_attr='intactMiscore') #  ,create_using=nx.MultiGraph()


    sys.stdout.write('\rDefining source and target nodes... ✔')

    print()
    sys.stdout.write('\rCreating layout...')
    plt.figure(figsize=(100, 100))
    plt.tight_layout()
    pos = nx.kamada_kawai_layout(G)

    sys.stdout.flush()
    sys.stdout.write('\rCreating layout... ✔')
    print()



    if labels == 'all' or labels == 'protein_interactions':
        sys.stdout.write('Creating labels...')
        sys.stdout.flush()
        sys.stdout.write('\rCreating labels... ✔')
        sys.stdout.flush()
        print()
        sys.stdout.write('Constructing network visualization...')
        nx.draw(G, font_color='red', node_color='lightblue', node_size=300, pos=pos, with_labels=True)
        sys.stdout.flush()
        sys.stdout.write('\rConstructing network visualization... ✔')
        sys.stdout.flush()
        print()


    else:
        sys.stdout.write('Creating network visualization...')
        sys.stdout.flush()
        nx.draw(G, font_color='red', node_color='lightblue', node_size=300, pos=pos, with_labels=False)
        sys.stdout.flush()
        sys.stdout.write('\rCreating network visualization... ✔')
        sys.stdout.flush()
        print()
num_nodes = G.number_of_nodes()


# ---------------------------------------------------------------------------
# Summary statistics calculation |
# -------------------------------+
print()
print('Calculating connectivity statistics for summary file:')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()

# Begin loading spinner animation
calculation_wait_animation.start()


components = nx.connected_components(G)
largest_component = max(components, key=len)

subgraph = nx.subgraph(G, largest_component)
diameter = nx.diameter(subgraph)

transitivity = nx.transitivity(G)

degree_dict = dict(nx.degree(G,G.nodes))
nx.set_node_attributes(G,degree_dict,'degree')
sorted_degree = sorted(degree_dict.items(),key=itemgetter(1),reverse=True)

betweenness_dict = nx.betweenness_centrality(G)
eigenvector_dict = nx.eigenvector_centrality(G,max_iter=1000)

nx.set_node_attributes(G,betweenness_dict,'betweenness')
nx.set_node_attributes(G,eigenvector_dict,'eigenvector')

sorted_betweenness = sorted(betweenness_dict.items(), key=itemgetter(1),reverse=True)

communities = nx.community.greedy_modularity_communities(G)
node_connectivity = nx.average_node_connectivity(G,flow_func=shortest_augmenting_path)
modularity_dict = {} # Create a blank dictionary
for i,c in enumerate(communities): # Loop through the list of communities, keeping track of the number for the community
    for name in c: # Loop through each person in a community
        modularity_dict[name] = i # Create an entry in the dictionary for the person, where the value is which group they belong to.

# Now you can add modularity information like we did the other metrics
nx.set_node_attributes(G, modularity_dict, 'modularity')



summary_title = """                                            
                                            __                               
   /| |      /                   /         /                                 
  ( | | ___ (___       ___  ___ (         (___       _ _  _ _  ___  ___      
  | | )|___)|    |   )|   )|   )|___)         )|   )| | )| | )|   )|   )\   )
  | |/ |__  |__  |/\/ |__/ |    | \        __/ |__/ |  / |  / |__/||     \_/ 
                                                                          /  
                                                 ₲Ɇ₦Ø₱ⱧɆ₦Ø v5.5          /
------------------------------------------------------------------------------
"""
summary_file = open(output.split('.')[0]+'_NETWORK_SUMMARY.txt', 'w')
print(summary_title, file = summary_file)
print('', file= summary_file)
print('Calculation date/time:  ',datetime.now(), file= summary_file)
print('', file= summary_file)
print('Total Nodes:  '+str(G.number_of_nodes()), file= summary_file)
print('Total Edges:  '+str(G.number_of_edges()), file= summary_file)
print('', file= summary_file)
print('                 Average node connectivity:  ', node_connectivity, file = summary_file)
print('     Network diameter of largest component:  ', diameter, file = summary_file)
print('                              Transitivity:  ', transitivity, file = summary_file)
print('', file= summary_file)
print('  ____        _             _      ', file = summary_file)
print('   L|op  20  [|\|odes  by  [|)egree ', file = summary_file)
print('   ---------------------------------', file = summary_file)
print('', file= summary_file)
print('Node  |  Degree', file= summary_file)
print('', file= summary_file)
sys.stdout.flush()
for d in sorted_degree[:20]:
    print(str(d[0])+':   '+str(d[1]), file= summary_file)

print('', file= summary_file)
print('', file= summary_file)
print('  ____        _             _              _ ', file = summary_file)
print('   L|op  20  [|\|odes  by  [|}etweenness  ((entrality', file = summary_file)
print('   --------------------------------------------------', file = summary_file)
print('', file= summary_file)
print('Node  |  Betweenness Centrality', file= summary_file)
print('', file= summary_file)
sys.stdout.flush()
for b in sorted_betweenness[:20]:
    print(str(b[0])+':   '+str(b[1]), file= summary_file)
print('', file= summary_file)

print("   _ _             _      ", file= summary_file)
print("  //\/\odularity  ((lasses (with top 10 nodes listed by Eigenvector centrality)", file= summary_file)
print(' -------------------------------------------------------------------------------', file = summary_file)
print('', file = summary_file)
print('Node | Eigenvector Centrality', file = summary_file)
print('', file = summary_file)
for i in range(len(communities)):
    # First get a list of just the nodes in that class
    class0 = [n for n in G.nodes() if G.nodes[n]['modularity'] == i]

    # Then create a dictionary of the eigenvector centralities of those nodes
    class0_eigenvector = {n: G.nodes[n]['eigenvector'] for n in class0}

    # Then sort that dictionary and print the first 5 results
    class0_sorted_by_eigenvector = sorted(class0_eigenvector.items(), key=itemgetter(1), reverse=True)
    print('', file=summary_file)
    print('CLASS '+str(i)+':', file = summary_file)
    print('-------', file=summary_file)
    for node in class0_sorted_by_eigenvector[:10]:
        print(str(node[0])+":   "+str(node[1]), file= summary_file)

# End loading bar
done = True
time.sleep(2)

sys.stdout.flush()


# ---------------------------------------------------------------------------
# Save output to files |
# ---------------------+
print()
print()
print()
print('Saving your output...')
print()
# Name the graph output file based on the input argument for the file name.
# Append '.png' to the filename and save the figure as that filename.
graph_output_name = output.split('.')[0]
graph_output_name += '.png'
plt.savefig(graph_output_name)

gexf_output_name = output.split('.')[0]
gexf_output_name += '.gexf'
nx.write_gexf(G, gexf_output_name)

# ---------------------------------------------------------------------------
# Print logo and output message |
# ------------------------------+

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


"""
print('', file = summary_file)
print('', file = summary_file)
print('Classes:', file = summary_file)
for i,c in enumerate(communities): # Loop through the list of communities
    if len(c) > 10: # Filter out modularity classes with 10 or fewer nodes
        print('Class '+str(i)+':', list(c), file= summary_file) # Print out the classes and their members
"""
summary_file.close()

if mode == 'geno':
    print(spiderweb_ascii_art)
if mode == 'pheno':
    print(spiderweb_ascii_art2)
print()
print()
print('     ...image saved as \"',output,'\" with ', num_nodes,' total nodes.')
print()
print('     Summary file saved as \"'+output.split('.')[0]+'_NETWORK_SUMMARY.txt\"')
print()
print('     Graph exchange XML file saved as \"'+gexf_output_name+'  *:･ﾟ✧')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print()
sys.stdout.flush()