#                              d8,        d8b
#                             `8P         ?88                                        d8P
#                                          88b                                    d888888P |
#   d8888b  88bd88b   88bd88b  88b d8888b  888888b   88bd8b,d88b  d8888b  88bd88b   ?88'   |
#  d8b_,dP  88P' ?8b  88P'  `  88Pd8P' `P  88P `?8b  88P'`?8P'?8bd8b_,dP  88P' ?8b  88P    |
#  88b     d88   88P d88      d88 88b     d88   88P d88  d88  88P88b     d88   88P  88b    |
#  `?888P'd88'   88bd88'     d88' `?888P'd88'   88bd88' d88'  88b`?888P'd88'   88b  `?8b   |
#                                                                                          |
#  ----------------------------------------------------------------------------------------+
#                                                        ₲Ɇ₦Ø₱ⱧɆ₦Ø  6.1
#
#                                                        Created: 2022-03-15
#
#  -----------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import argparse
import json
import requests
import ast
import matplotlib.pyplot as plt
from io import BytesIO
import itertools
import threading
import time
import sys
import os
import csv

logo = """
                              d8,        d8b
                             `8P         ?88                                        d8P
                                          88b                                    d888888P |
   d8888b  88bd88b   88bd88b  88b d8888b  888888b   88bd8b,d88b  d8888b  88bd88b   ?88'   |
  d8b_,dP  88P' ?8b  88P'  `  88Pd8P' `P  88P `?8b  88P'`?8P'?8bd8b_,dP  88P' ?8b  88P    |
  88b     d88   88P d88      d88 88b     d88   88P d88  d88  88P88b     d88   88P  88b    |
  `?888P'd88'   88bd88'     d88' `?888P'd88'   88bd88' d88'  88b`?888P'd88'   88b  `?8b   |
                                                                                          |
  ----------------------------------------------------------------------------------------+
                                                        ₲Ɇ₦Ø₱ⱧɆ₦Ø  v6.1
"""
print(logo)

# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Perform enrichment analysis on your protein interactors edge list! * ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_FILENAME.csv>  (protein interactors table from interactors.py \"geno\" mode)')
parser.add_argument('-m', '--mode', type=str, help='\"omim\" (table.py output) or \"interactors\" (interactors.py output)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME.csv>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
mode = args.mode


sys.stdout.flush()
print()
print('Converting gene IDs to Entrez:')
print('------------------------------')
print()
toppGene_command = ''

if mode == 'omim':
    gpn = pd.read_csv(args.input)
    gene_ids_list = []
    gene_ids_list_unique = []
    for i in range(len(gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols'].reset_index())):
        gene_ids_list += [gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols'].reset_index()['Node_name'][i]]

    for i in gene_ids_list:
        if i in gene_ids_list_unique:
            pass
        else:
            gene_ids_list_unique += [i]

    toppGene_command = 'curl -H \'Content-Type: text/json\' -d \'{\"Symbols\":['

    for i in range(len(gene_ids_list_unique)):
        toppGene_command += '\"'
        toppGene_command += gene_ids_list_unique[i]
        if i < len(gene_ids_list_unique)-1:
            toppGene_command += '\",'
        else:
            toppGene_command += '\"'

    toppGene_command += ']}\' https://toppgene.cchmc.org/API/lookup > '
    toppGene_command += 'id_conversion.json'

        #r = requests.post(toppGene_url)
        #response = r.json()
        #dictionary.update(response)

sys.stdout.flush()

if mode == 'interactors':
    interactors_df = pd.read_csv(input)

    interactors_list_unique = []

    interactors_list_A = list(interactors_df['moleculeA'])
    interactors_list_B = list(interactors_df['moleculeB'])



    interactors_list_total = interactors_list_A + interactors_list_B

    for i in interactors_list_total:
        if i in interactors_list_unique:
            pass
        else:
            interactors_list_unique += [i]

    toppGene_command = 'curl -H \'Content-Type: text/json\' -d \'{\"Symbols\":['

    for i in range(len(interactors_list_unique)):
        toppGene_command += '\"'
        toppGene_command += interactors_list_unique[i]
        if i < len(interactors_list_unique)-1:
            toppGene_command += '\",'
        else:
            toppGene_command += '\"'

    toppGene_command += ']}\' https://toppgene.cchmc.org/API/lookup > '
    toppGene_command += 'id_conversion.json'

        #r = requests.post(toppGene_url)
        #response = r.json()
        #dictionary.update(response)

os.system(toppGene_command)
sys.stdout.flush()
print()
data = pd.read_json('id_conversion.json')

entrez_list = []
for i in range(len(data['Genes'])):
    entrez_list += [data['Genes'][i]['Entrez']]

sys.stdout.flush()

print()
print('Performing enrichment analysis on Entrez IDs:')
print('---------------------------------------------')
print()

final_df = pd.DataFrame()

chunked_list = []
chunk_size = 1000

for i in range(0, len(entrez_list), chunk_size):
    chunked_list.append(entrez_list[i:i+chunk_size])


#  ** ** ADD LIMIT AND CHUNK LOOP HERE ** **
for j in range(len(chunked_list)):
    toppGene_command2 = 'curl -H \'Content-Type: text/json\' -d \'{\"Genes\":['
    for i in range(len(chunked_list[j])):
        toppGene_command2 += str(chunked_list[j][i])
        if i < len(chunked_list[j])-1:
            toppGene_command2 += ','
        else:
            pass

    toppGene_command2 += ']}\' https://toppgene.cchmc.org/API/enrich > ToppGene_response.json'
    sys.stdout.flush()
    os.system(toppGene_command2)
    sys.stdout.flush()
    data2 = pd.read_json('ToppGene_response.json')

    data3 = pd.json_normalize(data2['Annotations'])

    chunk_df = data3[data3['QValueFDRBH'] < 0.000001].sort_values(by='QValueFDRBH')

    final_df = pd.concat([final_df,chunk_df])

sys.stdout.flush()

deletion_command = 'rm ToppGene_response.json id_conversion.json'
os.system(deletion_command)
sys.stdout.flush()

if len(final_df) == 0:
    print('-------------------------------------------------------------------------------------')
    print()
    print('          Sorry, your enrichment analysis did not produce results meeting the')
    print('                     minimum threshold for statistical significance.')
    print()
    print('-------------------------------------------------------------------------------------')

else:
    print()
    print('-------------------------------------------------------------------------------------')
    print('-------------------------------------------------------------------------------------')
    print()
    print('Your output table:')

    final_df = final_df.reset_index(drop=True)
    final_df.to_csv(output)

    print(final_df)
    print()



    print('     * Your enrichment analysis was saved as \"'+output+'\"')
    print('     * Your table was filtered and sorted in ascending order by B&H/FDR Q Value < 10e-6')
    print()
    print('-------------------------------------------------------------------------------------')
    print()
    print('                        Thank you for using ₲Ɇ₦Ø₱ⱧɆ₦Ø    ❤')
    print()
    print('-------------------------------------------------------------------------------------')
    print()


"""
final_output = pd.DataFrame()
for i in range(len(data2)):
    tempdf =  pd.DataFrame.from_dict(data2[i])
    final_output = pd.concat([final_output,tempdf], axis=0, ignore_index=True)
"""