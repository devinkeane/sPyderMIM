#!/usr/bin/env python
# coding: utf-8

#     ___|                       _ \   |
#    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
#    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
#   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
#   ______________________________________________________________|
#                                        (ツ)_/¯   - * B E T A * -
#  (c) 2022-01-27 by Devin Keane
#  rev 1.1 --> 2022-02-02
#  Feltus Lab
#  Department of Genetics and Biochemistry, Clemson University

# ---------------------------------------------------------------------------
# Import libraries
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import networkx as nx
# ---------------------------------------------------------------------------
# Parse command line input and options
parser = argparse.ArgumentParser(description='Calculate the volume of a cylinder')
parser.add_argument('-i', '--input', type=str, help='MIM (OMIM) reference number list as .txt')
parser.add_argument('-a', '--apikey', type=str, help='MIM reference number (OMIM)')
parser.add_argument('-o', '--output', type=str, help='MIM reference number (OMIM)')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
# ---------------------------------------------------------------------------
# Populate a new dataframe with the input (.txt file of MIM numbers)
mimdf = pd.read_csv(input, header=None)

# Begin building the OMIM API request URL
url = 'https://api.omim.org/api/entry?mimNumber='

# Append the list of MIM numbers (now in "mimdf") onto URL
for line in range(len(mimdf)):
    url += mimdf[0][line].astype(str)
    url += ','

# Append clinical synopsis API request, JSON format option onto URL
url += '&include=clinicalSynopsis&format=json&apiKey='

# Append API key from parsed command line arguments onto URL
url += args.apikey

url_geneMap = url.replace('clinicalSynopsis','geneMap')

# ---------------------------------------------------------------------------
# Load data into new data frame using pandas .read_json() module
df = pd.read_json(url)

# flatten data
df2 = pd.json_normalize(df['omim'][0])
# - - - - - - - - - -
# load geneMap data using Python JSON module
df_geneMap = pd.read_json(url_geneMap)

# flatten data
df_geneMap2 = pd.json_normalize(df_geneMap['omim'][0])
df_geneMap2['entry.phenotypeMapList'] = pd.json_normalize(df_geneMap2['entry.phenotypeMapList'])

# drop molecular basis
# (not so data friendly gene id column, we will use the nested geneMap list from the
# API request instead)
df2.drop(columns='entry.clinicalSynopsis.molecularBasis',inplace=True)

# Transpose the clinical data --> now each column represents a disease (MIM#)
df2_transposed = pd.DataFrame.transpose(df2)

# Transpose the geneMap data --> now each column represents a disease (MIM#)
df_geneMap2_transposed = pd.DataFrame.transpose(df_geneMap2)

# Create a new data frame that will be the output and call it "gpn" (genotype/phenotype network)
gpn = pd.DataFrame(columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number'])

# ---------------------------------------------------------------------------
# The following function calculates how many total elements are found for a given
# column of df2_transposed without having to flatten all nested lists:
# ---------------------------------------------------------------------------
def count_elements(column):
    # Initialize the count
    count = 0

    # df2_transposed is structured in such a way that we need to begin on index 5
    # in order to access the phenotypic data that we will populate gpn with
    i = 5

    # For each row in df2_transposed, starting at row 5:
    for row in range(len(df2_transposed) - 5):

        # If the current cell is null/float/integer, nothing happens to the count
        if pd.isnull(df2_transposed[column].iloc[i]) | pd.api.types.is_float(
                df2_transposed[column].iloc[i]) | pd.api.types.is_integer(df2_transposed[column].iloc[i]):
            i += 1

        # If the cell is not null but no '\n' is found, there is only one element in the cell
        # Thus, only add one to the count
        elif df2_transposed[column].iloc[i].count('\n') == 0:
            count += 1
            i += 1

        # Otherwise, the amount of elemnts should be equal to n + 1
        elif df2_transposed[column].iloc[i].count('\n') > 0:
            count += df2_transposed[column].iloc[i].count('\n') + 1
            i += 1

    return count
# ---------------------------------------------------------------------------------------------
# Populating a fresh genotype/phenotype dataframe with data from df2_transposed (OMIM API data)
# ---------------------------------------------------------------------------------------------
# We need to find the total amount of items in df2_transposed in order to build
# a new data frame that will ultimately be the genotype/phenotype table, which
# we will call "gpn".

# Count the total amount of elements across by calling count_elements() on each column
total_elements = 0
for i in range(len(df2_transposed.columns)):
    total_elements += count_elements(i)

# Now we can prepare gpn as a data frame of n rows, where n = total_elements
for i in range(total_elements):
    dfa = pd.DataFrame([[np.nan] * len(gpn.columns)], columns=gpn.columns)
    gpn = dfa.append(gpn, ignore_index=True)

# Create the variable k, which will represent the iteration of the index of gpn
# as we populate it with data from df2_transposed.
k = 0

# For each column in d2_transposed (each column is a disease)...
for j in range(len(df2_transposed.columns)):

    # take each element in that column and add it to the 'Superphenotype' column
    # of gpn at row k, but only if something is there (not null)
    for i in range(count_elements(j)):
        if pd.notnull(df2_transposed[j][3]):
            gpn['Superphenotype'][k] = df2_transposed[j][3].split('; ')[-1]
            k += 1

# Reset indices:  go back to the top row of gpn and row at index 5 on df2_transposed
k = 0
i = 5

# for every row of every column in df2_transposed, starting on row = 5...
for j in range(len(df2_transposed.columns)):
    for row in range(len(df2_transposed)):
        if row > 4:

            # If the current cell is null/float/integer, skip it
            if pd.isnull(df2_transposed[j][row]) | pd.api.types.is_float(df2_transposed[j][row]) | pd.api.types.is_integer(df2_transposed[j][row]):
                 i += 1

            # If the cell is not null but no '\n' is found, there is only one element in the cell.
            # Add this element to the 'Node_name' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') == 0:
                gpn['Node_name'][k] = df2_transposed[j][row]
                gpn['Node_type'][k] = df2_transposed.index[row]
                gpn['MIM_number'][k] = df_geneMap2_transposed[j]['entry.mimNumber']
                k += 1
                i += 1

            # Otherwise, the amount of elemnts should be equal to n + 1,
            # so add each of these elements to the 'Node_name' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') > 0:
                for element in range(len(df2_transposed[j][row].split('\n'))):
                    gpn['Node_name'][k] = df2_transposed[j][row].split('\n')[element]
                    gpn['Node_type'][k] = df2_transposed.index[row]
                    gpn['MIM_number'][k] = df_geneMap2_transposed[j]['entry.mimNumber']
                    k += 1
                i += 1
# ---------------------------------------------------------------------------
# Now, we need to clean up the Node_name column data in gpn by cutting the
# off the indentifiers and other data that follows it.

# So, for each row in gpn...
for i in range(len(gpn)):
    # if it is a list, ignore it
    if isinstance(gpn['Node_name'][i], list):
        pass
    # Otherwise, if it is a string and '{' can be found in it, cut that off and everything
    # that comes after it as well.
    else:
        if isinstance(gpn['Node_name'][i], str):
            if gpn['Node_name'][i].count('{') > 0:
                gpn['Node_name'][i] = gpn['Node_name'][i].split(' {', 1)[0]
# ---------------------------------------------------------------------------
# GRAPHING THE DATA

# Create a NetworkX object called "G" where 'Superphenotype' is the source node
# and 'Node_name' is the target node.
G = nx.from_pandas_edgelist(gpn,source = 'Superphenotype', target = 'Node_name')

# Draw a graph with G, setting node color to green and size to 7
nx.draw(G, node_size=7, node_color='green')

# Name the graph output file based on the input argument for the file name.
# Append '.png' to the filename and save the figure as that filename.
graph_output_name = output.split('.')[0]
graph_output_name += '.png'
plt.savefig(graph_output_name)

# Save the gpn as a csv using the same filename, but with extension '.csv'
gpn.to_csv(output)
# ---------------------------------------------------------------------------
# Print logo and output message
# ------------------------------
logo = """
     ___|                       _ \   |
    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
   ______________________________________________________________|
                           ⌒ *: ﾟ･✧* ･ﾟ✧ - * [1.1] 🅱 🅴 🆃 🅰 * -
                    (ツ)_/¯
"""
print()
print('Thank you for using...')
print(logo)
print('                 ...Your network table was saved as \"',output,'\" with ',len(gpn),' total rows.')
print('                    Your network summary graph was saved as \"',graph_output_name,'\".')
print()

# End of program