#                               ,--.
#   ,---. ,--.--. ,--,--. ,---. |  ,---.      ,---.,--. ,--.
#  | .-. ||  .--'' ,-.  || .-. ||  .-.  |    | .-. |\  '  /
#  ' '-' '|  |   \ '-'  || '-' '|  | |  |.--.| '-' ' \   '
#  .`-  / `--'    `--`--'|  |-' `--' `--''--'|  |-'.-'  /
#  `---'                 `--'                `--'  `---'
#                                [  s P y d e r M I M  ]
#
# last rev: 2022-10-05
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
parser.add_argument('-m', '--mode', type=str, required=True , help='\"interactors\" (find protein interactor overlap), \"omim_genes\" (find OMIM genes shared between OMIM diseases), \"omim_genes_interactions\" (find interactions between OMIM genes), or  \"omim_features\" (find OMIM clinical feature overlap)')
parser.add_argument('-i', '--input', type=str, required=True ,help='<INPUT_FILENAME.csv>  (Input table)')
parser.add_argument('-g', '--gene_list', required=False, type=str, default='none', help='<GENES_LIST.csv>  OMIM Genes List Output (only for \'omim_genes_interactions\' mode)')
parser.add_argument('-l', '--labels', required=False, type=str, default='none', help='Labels --> arguments: \"all\" or \"none\"')
parser.add_argument('-o', '--output', type=str, required=True , help='<OUTPUT_NAME> (without file extension)')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
labels = args.labels
gene_list = args.gene_list
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

# ---------------------------------------------------------------------------
# TITLE SCREEN / LOGO |
# --------------------+

logo = """

                    ____            __          __  _________  ___   |
              _____/ __ \__  ______/ /__  _____/  |/  /  _/  |/  /   |
O---o        / ___/ /_/ / / / / __  / _ \/ ___/ /|_/ // // /|_/ /    |
O---o       (__  ) ____/ /_/ / /_/ /  __/ /  / /  / // // /  / /     |
 O-o       /____/_/    \__, /\__,_/\___/_/  /_/  /_/___/_/  /_/      |
  O                   /____/                                         |
 o-O    _____________________________________________________________|----------+
o---O   High Performance Computing Genomic Network Analysis    |  Version 7.4.0 |   ✧ - ･ﾟ*
O---o   :". /  /  /                          +------------------------------------------+ 
 O-o    :.-". /  /                           |  (c) 2022-01-27 Devin Keane              |
  O     : _.-". /                            |  Feltus Lab                              |
 o-O    :"  _.-".                            |  Department of Genetics and Biochemistry |
o---O   :-""     ".                          |  Clemson University                      | 
O---o   :                          .'✧       |                                          |
        :                                    |  Last rev: 2023-05-06                    |
      ^.-.^                                  +------------------------------------------+
     '^\+/^`
     '/`"'\`

"""
print(logo)

# ---------------------------------------------------------------------------
# GRAPHING THE DATA |
# ------------------+

G = nx.Graph()

if mode == 'omim_genes':
    gpn = pd.read_csv(input)
    print('Processing your input table:')
    print()
    print(gpn.drop(columns='Unnamed: 0'))
    print()

    gpn = gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols']
    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'Superphenotype' is the source node
    # and 'Node_name' is the target node.
    G = nx.from_pandas_edgelist(gpn,source = 'Superphenotype', target = 'Node_name')
    done = False
    while done == False:
        nodeList = list(G.nodes)
        count = 0
        for node in nodeList:
            if G.degree(node) <= 1:
                if node not in list(gpn['Superphenotype']):
                    G.remove_node(node)
                    count += 1
        if count == 0:
            done = True


    nodeList = list(G.nodes)
    for node in nodeList:
        if G.degree(node) == 0:
            G.remove_node(node)

    nodeList = list(G.nodes)
    for node in nodeList:
        try:
            if G[node][node] != None:
                if G.degree(node) <= 2:
                    G.remove_node(node)
        except Exception:
            pass

    sys.stdout.flush()

    sys.stdout.flush()

    sys.stdout.write('\rDefining source and target nodes... ✔')
    print()

    sys.stdout.write('Creating layout...')
    # Draw a graph with G using a color map that distinguishes between genes and phenotypes
    plt.figure(figsize=((G.number_of_nodes()/2)+(5),(G.number_of_nodes()/2)+(5)))
    sys.stdout.flush()






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

    nx.draw(G, node_color='green', pos=pos,with_labels=False)


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

        bbox = dict(fc="blue", ec="black", boxstyle="square")



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
if mode == 'omim_features':
    gpn = pd.read_csv(input)
    print('Processing your input table:')
    print()
    print(gpn.drop(columns='Unnamed: 0'))
    print()

    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'Superphenotype' is the source node
    # and 'Node_name' is the target node.
    G = nx.from_pandas_edgelist(gpn, source='Superphenotype', target='Node_name')

    done = False
    while done == False:
        nodeList = list(G.nodes)
        count = 0
        for node in nodeList:
            if G.degree(node) <= 1:
                if node not in gpn['Superphenotype']:
                    G.remove_node(node)
                    count += 1
                count += 1
        if count == 0:
            done = True

    nodeList = list(G.nodes)
    for node in nodeList:
        if G.degree(node) == 0:
            G.remove_node(node)

    nodeList = list(G.nodes)
    for node in nodeList:
        try:
            if G[node][node] != None:
                if G.degree(node) <= 2:
                    G.remove_node(node)
        except Exception:
            pass

    sys.stdout.flush()
    sys.stdout.write('\rDefining source and target nodes... ✔')
    print()

    sys.stdout.write('Creating layout...')
    # Draw a graph with G using a color map that distinguishes between genes and phenotypes
    #plt.figure(figsize=(100, 100))
    sys.stdout.flush()
    plt.figure(figsize=((G.number_of_nodes()/2)+(5),(G.number_of_nodes()/2)+(5)))
    pos = nx.kamada_kawai_layout(G)
    # pos= nx.spring_layout(G)
    sys.stdout.flush()
    sys.stdout.write('\rCreating layout... ✔')
    print()

    sys.stdout.write('Constructing network visualization...')

    nx.draw(G, node_color='green', pos=pos, with_labels=False)

    # color_map = []
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

        bbox = dict(fc="blue", ec="black", boxstyle="square")

        nx.draw_networkx_labels(G, pos, labels=superphenotype_labels, font_color='white',
                                font_family='copperplate', bbox=bbox)
        nx.draw_networkx_labels(G, pos, labels=Node_name_labels, font_color='black',
                                font_family='copperplate')

        sys.stdout.write('\rCreating labels... ✔')
        sys.stdout.flush()
        print()

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


elif mode == 'interactors':
    protein_df = pd.read_csv(input)
    print('Processing your input table:')
    print()
    print(protein_df.drop(columns='Unnamed: 0'))

    print()
    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'moleculeA' is the source node
    # and 'moleculeB' is the target node.
    G = nx.from_pandas_edgelist(protein_df, source='moleculeA', target='moleculeB',edge_attr='intactMiscore') #  ,create_using=nx.MultiGraph()

    done = False
    while done == False:
        nodeList = list(G.nodes)
        count = 0
        for node in nodeList:
            if G.degree(node) <= 1:
                G.remove_node(node)
                count += 1
        if count == 0:
            done = True

    nodeList = list(G.nodes)
    for node in nodeList:
        try:
            if G[node][node] != None:
                if G.degree(node) <= 2:
                    G.remove_node(node)
        except Exception:
            pass
    # The following code may be used in the future if a subplot is desired that would remove all
    # nodes of degree n=1.
    """
    to_be_removed = [x for x in G.nodes() if G.degree(x) <= 1]

    for x in to_be_removed:
        G.remove_node(x)
    """

    sys.stdout.write('\rDefining source and target nodes... ✔')

    print()
    sys.stdout.write('\rCreating layout...')
    #plt.figure(figsize=(100, 100))
    plt.figure(figsize=((G.number_of_nodes() / 2) + (5), (G.number_of_nodes() / 2) + (5)))

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
        nx.draw(G, font_color='red', node_color='lightblue',pos=pos, with_labels=True)
        sys.stdout.flush()
        sys.stdout.write('\rConstructing network visualization... ✔')
        sys.stdout.flush()
        print()

    else:
        sys.stdout.write('Creating network visualization...')
        sys.stdout.flush()
        nx.draw(G, font_color='red', node_color='lightblue', pos=pos, with_labels=False)
        sys.stdout.flush()
        sys.stdout.write('\rCreating network visualization... ✔')
        sys.stdout.flush()
        print()


if mode == 'omim_genes_interactions':
    protein_df = pd.read_csv(input)
    gene_list_df = pd.read_csv(gene_list)

    print('Processing your input table:')
    print()
    print(protein_df.drop(columns='Unnamed: 0'))

    protein_df.drop(columns='Unnamed: 0', axis=1, inplace=True)
    genes_df = pd.read_csv(gene_list)

    genes_list = []
    genes_list = genes_df['Node_name'][genes_df.Node_type == 'phenotypeMap.approvedGeneSymbols'].tolist()

    for i in range(len(protein_df.index)):
        if protein_df['moleculeA'][i] in genes_list:
            if protein_df['moleculeB'][i] in genes_list:
                pass
            else:
                protein_df.drop(i, inplace=True)
        else:
            protein_df.drop(i, inplace=True)

    protein_df.reset_index(drop=True,inplace=True)



    print()
    sys.stdout.write('Defining source and target nodes...')

    # Create a NetworkX object called "G" where 'moleculeA' is the source node
    # and 'moleculeB' is the target node.
    G = nx.from_pandas_edgelist(protein_df, source='moleculeA', target='moleculeB')  # ,create_using=nx.MultiGraph()
    nodeList = list(G.nodes)
    for node in nodeList:
        try:
            if G[node][node] != None:
                if G.degree(node) <= 2:
                    G.remove_node(node)
        except Exception:
            pass

    """
    for i in genes_list:
        if i not in G.nodes:
            G.add_node(i,pos=nx.kamada_kawai_layout(i))
    """

    sys.stdout.write('\rDefining source and target nodes... ✔')

    print()
    sys.stdout.write('\rCreating layout...')
    #plt.figure(figsize=(25, 25))
    plt.figure(figsize=((G.number_of_nodes() / 2) + (5), (G.number_of_nodes() / 2) + (5)))

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
        nx.draw(G, font_color='red', node_color='lightblue', pos=pos, with_labels=True)
        sys.stdout.flush()
        sys.stdout.write('\rConstructing network visualization... ✔')
        sys.stdout.flush()
        print()

    else:
        sys.stdout.write('Creating network visualization...')
        sys.stdout.flush()
        nx.draw(G, font_color='red', node_color='lightblue', pos=pos, with_labels=False)
        sys.stdout.flush()
        sys.stdout.write('\rCreating network visualization... ✔')
        sys.stdout.flush()
        print()

# ---------------------------------------------------------------------------
# Summary statistics calculation |
# -------------------------------+
print()
print('Calculating connectivity statistics for summary file:')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()

if not any(G.adj.values()):
    summary_title = """                                            
                                                __                               
       /| |      /                   /         /                                 
      ( | | ___ (___       ___  ___ (         (___       _ _  _ _  ___  ___      
      | | )|___)|    |   )|   )|   )|___)         )|   )| | )| | )|   )|   )\   )
      | |/ |__  |__  |/\/ |__/ |    | \        __/ |__/ |  / |  / |__/||     \_/ 
                                                                              /  
                                                     𝕤ℙ𝕪𝕕𝕖𝕣𝕄𝕀𝕄 v7.4.0         /
    ------------------------------------------------------------------------------
    """
    summary_file = open(output + '_NETWORK_SUMMARY.txt', 'w')
    print(summary_title, file=summary_file)
    print('', file=summary_file)
    print('Calculation date/time:  ', datetime.now(), file=summary_file)
    print('', file=summary_file)
    print('Total Nodes:  0', file=summary_file)
    print('Total Edges:  0', file=summary_file)
    print('', file=summary_file)
    print('             Average degree connectivity:  0', file=summary_file)
    print('     Network diameter of largest component:  0', file=summary_file)
    print('                              Transitivity:  0', file=summary_file)
    print('', file=summary_file)
else:
    num_nodes = G.number_of_nodes()
    # Begin loading spinner animation
    #calculation_wait_animation.start()
    components = nx.connected_components(G)

    largest_component = max(components, key=len)

    subgraph = nx.subgraph(G, largest_component)
    diameter = nx.diameter(subgraph)

    transitivity = nx.transitivity(G)

    degree_dict = dict(nx.degree(G,G.nodes))
    nx.set_node_attributes(G,degree_dict,'degree')
    sorted_degree = sorted(degree_dict.items(),key=itemgetter(1),reverse=True)

    betweenness_dict = nx.betweenness_centrality(G)
    eigenvector_dict = nx.eigenvector_centrality(G,max_iter=10000)

    nx.set_node_attributes(G,betweenness_dict,'betweenness')
    nx.set_node_attributes(G,eigenvector_dict,'eigenvector')

    sorted_betweenness = sorted(betweenness_dict.items(), key=itemgetter(1),reverse=True)

    communities = nx.community.greedy_modularity_communities(G)
    node_connectivity = 2*(G.number_of_edges()/G.number_of_nodes())
    # node_connectivity = nx.average_node_connectivity(G,flow_func=shortest_augmenting_path)
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
                                                     𝕤ℙ𝕪𝕕𝕖𝕣𝕄𝕀𝕄 v7.3.1         /
    ------------------------------------------------------------------------------
    """
    summary_file = open(output+'_NETWORK_SUMMARY.txt', 'w')
    print(summary_title, file = summary_file)
    print('', file= summary_file)
    print('Calculation date/time:  ',datetime.now(), file= summary_file)
    print('', file= summary_file)
    print('Total Nodes:  '+str(G.number_of_nodes()), file= summary_file)
    print('Total Edges:  '+str(G.number_of_edges()), file= summary_file)
    print('', file= summary_file)
    print('             Average degree connectivity:  ', format(node_connectivity,'.3E'), file = summary_file)
    print('     Network diameter of largest component:  ', diameter, file = summary_file)
    print('                              Transitivity:  ', format(transitivity,'.3E'), file = summary_file)
    print('', file= summary_file)
    print('  ____        _             _      ', file = summary_file)
    print('   L|op  20  [|\|odes  by  [|)egree ', file = summary_file)
    print('   ---------------------------------', file = summary_file)
    print('', file= summary_file)
    print('{:>50} | {:>10}'.format('Node','Degree'), file= summary_file)
    print('', file= summary_file)
    sys.stdout.flush()
    for d in sorted_degree[:20]:
        print('{:>50}   {:>10}'.format(str(d[0]), str(d[1])), file= summary_file)

    print('', file= summary_file)
    print('', file= summary_file)
    print('  ____        _             _              _ ', file = summary_file)
    print('   L|op  20  [|\|odes  by  [|}etweenness  ((entrality', file = summary_file)
    print('   --------------------------------------------------', file = summary_file)
    print('', file= summary_file)
    print('{:>50} | {:>10}'.format('Node','Betweenness Centrality'), file= summary_file)
    print('', file= summary_file)
    sys.stdout.flush()
    for b in sorted_betweenness[:20]:
        print('{:>50}   {:>10}'.format(str(b[0]), format(b[1],'.3E')), file= summary_file)
    print('', file= summary_file)

    print("   _ _             _      ", file= summary_file)
    print("  //\/\odularity  ((lasses (with top 10 nodes listed by Eigenvector centrality)", file= summary_file)
    print(' -------------------------------------------------------------------------------', file = summary_file)
    print('', file = summary_file)
    print('{:>50} | {:>10}'.format('Node','Eigenvector Centrality'), file = summary_file)
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
            print('{:>50}   {:>10}'.format(str(node[0]),format(node[1],'.3E')), file= summary_file)
"""
# End loading bar
done = True
time.sleep(2)

sys.stdout.flush()
"""

# ---------------------------------------------------------------------------
# Save output to files |
# ---------------------+
print()
print()
print()
if not any(G.adj.values()):
    pass
else:
    print('Saving your output...')
    print()
    # Name the graph output file based on the input argument for the file name.
    # Append '.png' to the filename and save the figure as that filename.
    graph_output_name = output+'.png'

    plt.savefig(graph_output_name)

    graphml_output_name = output+'.graphml'
    nx.write_graphml(G, graphml_output_name)

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
             [  O M I M  ]                                       \  \/ | \/`. |___
                                                                  \ /`.|,'\  /
                                                                   Y.  |   \/
                                                                   | `.|_,'
                                                                   |
                                                                   |
                  Network spun.                                 __ |
                                                                __\|,-
                                                                ,-`=--.       
                                                                 /=8\\
"""

spiderweb_ascii_art2 = """
                                                                 \   ,'|`-.   /
                                                                  \,' _|_  ','
                                                                  /'.' | `,' \\
                                                              -._/_/_'.|,'_\__\_,-
                                                                 | | ,-*." |  |
                                                              ___|,+' /|\`.|  |
                        [  I n t e r a c t o r s  ]              \  \/ | \/`. |___
                                                                  \ /`.|,'\  /
                                                                   Y.  |   \/
                                                                   | `.|_,'
                                                                   |
                                                                   |
                  Network spun.                                 __ |
                                                                __\|,-
                                                                ,-`=--.       
                                                                 /=8\\
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

if mode == 'interactors':
    print(spiderweb_ascii_art2)
if mode == 'omim_genes' or mode == 'omim_features':
    print(spiderweb_ascii_art)
if not any(G.adj.values()):
    print('      ⚠ NO GRAPH WAS PRODUCED FOR THIS INPUT! ⚠')
    print()
    print('      No connected nodes were found.')
    print()
else:
    print()
    print()
    print('     ...image saved as \"'+output+'.png\" with ', num_nodes,' total nodes.')
    print()

    print('     Summary file saved as \"'+output+'_NETWORK_SUMMARY.txt\"')
    print()
    print('     GraphML file saved as \"'+graphml_output_name+'  *:･ﾟ✧')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print()
sys.stdout.flush()