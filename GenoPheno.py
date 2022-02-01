#!/usr/bin/env python
# coding: utf-8

#     ___|                       _ \   |
#    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
#    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
#   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
#   ______________________________________________________________|
#                                        (ツ)_/¯   - * B E T A * -
#  (c) 2022-01-27 by Devin Keane
#  Feltus Lab
#  Department of Genetics and Biochemistry, Clemson University

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import networkx as nx



# import json
# import csv

parser = argparse.ArgumentParser(description='Calculate the volume of a cylinder')
parser.add_argument('-i', '--input', type=str, help='MIM (OMIM) reference number list as .txt')
parser.add_argument('-a', '--apikey', type=str, help='MIM reference number (OMIM)')
parser.add_argument('-o', '--output', type=str, help='MIM reference number (OMIM)')
args = parser.parse_args()


input = args.input
output = args.output
mimdf = pd.read_csv(input, header=None)


url = 'https://api.omim.org/api/entry?mimNumber='
for line in range(len(mimdf)):
    url += mimdf[0][line].astype(str)
    url += ','

url += '&include=clinicalSynopsis&format=json&apiKey='
url += args.apikey

# consider adding this back to URL if needed in the future --> &include=geneMap

# url = 'https://api.omim.org/api/entry?mimNumber=615349,225400,130050,130010,614170,130070,225410,606408,616471,615539,614557,618000,225320,617821,130000,612350,617174,130080,601776,229200,130060,130020&include=geneMap&include=clinicalSynopsis&apiKey=JWd5fb0yR1a3TOUFrqE3UA&format=json'

# load data using Python JSON module
df = pd.read_json(url)
# flatten data
df2 = pd.json_normalize(df['omim'][0])
df2_transposed = pd.DataFrame.transpose(df2)

gpn = pd.DataFrame(columns=['Superphenotype', 'Node_name', 'Node_type'])

# The following function counts how many nested items are in a given column in df2_transposed


def count_elements(column):
    count = 0
    # We have to begin on the i where the phenotypic data starts
    i = 5

    for row in range(len(df2_transposed) - 5):

        # If the current cell is null, nothing happens to the count
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


total_nested = 0
k = 0

for i in range(len(df2_transposed.columns)):
    total_nested += count_elements(i)

for j in range(total_nested):
    dfa = pd.DataFrame([[np.nan] * len(gpn.columns)], columns=gpn.columns)
    gpn = dfa.append(gpn, ignore_index=True)


k = 0

for j in range(len(df2_transposed.columns)):
    for i in range(count_elements(j)):
        if pd.notnull(df2_transposed[j][3]):
            gpn['Superphenotype'][k] = df2_transposed[j][3].split('; ')[-1]
            k += 1

k = 0
i = 5

for j in range(len(df2_transposed.columns)):
    for row in range(len(df2_transposed)):
        if row > 4:
            # If the current cell is null, skip it
            if pd.isnull(df2_transposed[j][row]) | pd.api.types.is_float(df2_transposed[j][row]) | pd.api.types.is_integer(df2_transposed[j][row]):
                 i += 1
            # If the cell is not null but no '\n' is found, there is only one element in the cell
            elif df2_transposed[j][row].count('\n') == 0:
                gpn['Node_name'][k] = df2_transposed[j][row]
                gpn['Node_type'][k] = df2_transposed.index[row]
                k += 1
                i += 1
            # Otherwise, the amount of elemnts should be equal to n + 1
            elif df2_transposed[j][row].count('\n') > 0:
                for element in range(len(df2_transposed[j][row].split('\n'))):
                    gpn['Node_name'][k] = df2_transposed[j][row].split('\n')[element]
                    gpn['Node_type'][k] = df2_transposed.index[row]
                    k += 1
                i += 1
for i in range(len(gpn)):
    if isinstance(gpn['Node_name'][i], list):
        pass
    else:
        if isinstance(gpn['Node_name'][i], str):
            if gpn['Node_name'][i].count('{') > 0:
                gpn['Node_name'][i] = gpn['Node_name'][i].split(' {', 1)[0]

G = nx.from_pandas_edgelist(gpn,source = 'Superphenotype', target = 'Node_name')
nx.draw(G, node_size=7, node_color='green')

graph_output_name = output.split('.')[0]
graph_output_name += '.png'
plt.savefig(graph_output_name)
gpn.to_csv(output)
logo = """
     ___|                       _ \   |
    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
   ______________________________________________________________|

       (ツ)_/¯                            1.0 - * B E T A * -
"""
print()
print()
print('Thank you for using...')
print(logo)
print()
print()
print('                 ...Your network table was saved as \"',output,'\" with ',len(gpn),' total rows.')
print('                    Your network summary graph was saved as \"',graph_output_name,'\".')
print()
