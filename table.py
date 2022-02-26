#!/usr/bin/env python
# coding: utf-8

#     ___|                       _ \   |
#    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
#    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
#   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
#   ______________________________________________________________|
#                                      (ツ)_/¯  - * Version 3.0 * -
#  [ O m i m   T a b l e   M a k e r ]
#
# Last rev: 2022-02-26
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Import libraries
import numpy as np
import pandas as pd
import argparse

# ------------------------------------------------------------------------------------------------------
# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Build a genotype/phenotype network using an OMIM API key and a simple list of OMIM reference IDs! * ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_FILENAME.txt>  (list of MIM reference numbers, no headers, each MIM separated by a new line)')
parser.add_argument('-a', '--apikey', type=str, help='paste your API key obtained from OMIM ( https://www.omim.org/api )')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME.csv>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
# ------------------------------------------------------------------------------------------------------
# Print opening screen  |
# ----------------------+

logo = """
     ___|                       _ \   |
    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
   ______________________________________________________________|
                                      (ツ)_/¯  - * Version 3.0 * -
  [ O m i m   T a b l e   M a k e r ]
"""
print(logo)
# ------------------------------------------------------------------------------------------------------
# Populate a new dataframe with the input (.txt file of MIM numbers)
mimdf = pd.read_csv(input, header=None)
# ------------------------------------------------------------------------------------------------------
# Error --> MIM Exceeded Max |
# ---------------------------+
# If mim.txt contains over 20 mims, it cannot be run due to API request limits
if len(mimdf) > 20:
    print()
    print()
    print('Sorry, your query contains over 20 MIM numbers.   (ง ͠° ͟ʖ ͡°)')
    print('Please try again with a maximum of 20 MIM numbers')
    print()
    exit(-1)
# ------------------------------------------------------------------------------------------------------
# Inserting the MIM numbers into the API request URL and retrieving data |
# -----------------------------------------------------------------------+

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

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

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



# ------------------------------------------------------------------------------------------------------
# The following function calculates how many total elements are found for a  |
# given column of df2_transposed without having to flatten all nested lists: |
# ---------------------------------------------------------------------------+
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
# -----------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------+
# Populating a fresh genotype/phenotype dataframe with data from df2_transposed (OMIM API data) |
# ----------------------------------------------------------------------------------------------+

# Create a new data frame that will be the output and call it "gpn" (genotype/phenotype network)
gpn = pd.DataFrame(columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

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
    gpn = pd.concat([gpn,dfa], ignore_index=True)

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
            # Add this element to the 'Node_name_temp' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') == 0:
                gpn['Node_name_temp'][k] = df2_transposed[j][row]
                gpn['Node_type'][k] = df2_transposed.index[row]
                gpn['MIM_number'][k] = df_geneMap2_transposed[j]['entry.mimNumber']
                k += 1
                i += 1

            # Otherwise, the amount of elements should be equal to n + 1,
            # so add each of these elements to the 'Node_name_temp' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') > 0:
                for element in range(len(df2_transposed[j][row].split('\n'))):
                    gpn['Node_name_temp'][k] = df2_transposed[j][row].split('\n')[element].replace(';','')
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
    if isinstance(gpn['Node_name_temp'][i], list):
        pass
    # Otherwise, if it is a string and '{' can be found in it, cut that off and everything
    # that comes after it as well.
    else:
        if isinstance(gpn['Node_name_temp'][i], str):
            if gpn['Node_name_temp'][i].count('{') > 0:
                gpn['Node_name_temp'][i] = gpn['Node_name_temp'][i].split(' {', 1)[0].replace(';', '')
# ---------------------------------------------------------------------------
# Create parenthetical column
for i in range(len(gpn)):
    if gpn['Node_name_temp'][i].find('(') > 0:
        if gpn['Node_name_temp'][i].find('({') > 0:
            gpn['Node_name'][i] = gpn['Node_name_temp'][i]
        else:
            temp = gpn['Node_name_temp'][i].split(' (')[-1]
            gpn['Parenthetical'][i] = '('+temp

# Eliminate parentheticals, create permanent 'Node_name' column
for i in range(len(gpn)):
    if gpn['Node_name_temp'][i].find('(') > 0:
        if gpn['Node_name_temp'][i].find('({') > 0:
            if gpn['Node_name_temp'][i].find(') ') == 0:
                gpn['Node_name'][i] = gpn['Node_name_temp'][i]
        else:
            gpn['Node_name'][i] = gpn['Node_name_temp'][i].split(' (')[0]
    else:
        gpn['Node_name'][i] = gpn['Node_name_temp'][i]
gpn.drop(columns='Node_name_temp',inplace=True)

ensembl_df = pd.DataFrame(index=range(len(mimdf)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])


gpn2 = pd.DataFrame(index=range(len(mimdf)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

# Parsing ENSMBL IDs from df_geneMap2_transposed and other data from df2_transposed to write into gpn2
for i in range(len(mimdf)):
    gpn2['Node_name'][i] = df_geneMap2_transposed[i]['entry.phenotypeMapList']['phenotypeMap.ensemblIDs'].split(',')[0]
    gpn2['Node_type'][i] = 'phenotypeMap.ensemblIDs'
    gpn2['Superphenotype'][i] = df2_transposed[i][3].split('; ')[-1]
    gpn2['MIM_number'][i] = df_geneMap2_transposed[i]['entry.mimNumber']

gpn3 = pd.DataFrame(index=range(len(mimdf)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

for i in range(len(mimdf)):
    gpn3['Node_name'][i] = df_geneMap2_transposed[i]['entry.phenotypeMapList']['phenotypeMap.approvedGeneSymbols']
    gpn3['Node_type'][i] = 'phenotypeMap.approvedGeneSymbols'
    gpn3['Superphenotype'][i] = df2_transposed[i][3].split('; ')[-1]
    gpn3['MIM_number'][i] = df_geneMap2_transposed[i]['entry.mimNumber']

# Append gpn3 and gpn2 to gpn
gpn = pd.concat([gpn,gpn2], ignore_index=True)
gpn = pd.concat([gpn,gpn3], ignore_index=True)

ensembl_ids = (gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs'])['Node_name'].reset_index(drop=True)


# Save the gpn as a csv using the same filename, but with extension '.csv'
gpn.to_csv(output)
#df2_transposed.to_csv('testing.csv')

# Print output message
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print('  ...Your network table was saved as \"',output,'\" with ',len(gpn),' total rows using ',len(mimdf),' MIMs.')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()

# End of program

"""
#  ---------------------------------------------------------------------------------------
#  ---------------------------------------------------------------------------------------
#   **  The portion of code below has been removed but will remain here commented out   **
#   **  since it might very well be useful in future development.  -DK :)  2022-02-26   **
#  ---------------------------------------------------------------------------------------
#  ---------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# Intact Material | ** This section is being transferred to a separate program of its own **
#-----------------+


intact_url = ''
genes_and_overlap = []
for i in range(len(ensembl_ids)):
    intact_url = 'https://rest.ensembl.org/overlap/id/'
    intact_url += ensembl_ids[i]
    intact_url += '?content-type=application/json;feature=gene'
    print(intact_url)
    temp = pd.read_json(intact_url)
    for j in range(len(temp.external_name)):
        if pd.isnull(temp.external_name[j]):
            pass
        else:
            genes_and_overlap.append(temp.external_name[j])



# For testing: (remove later)
# for i in genes_and_overlap:
#    print(i)


genes_and_overlap_node_name = []
genes_and_overlap_node_type = []
genes_and_overlap_superphenotype = []
genes_and_overlap_mim_number = []

intact_url = ''

for i in range(len(ensembl_ids)):
    intact_url = 'https://rest.ensembl.org/overlap/id/'
    intact_url += ensembl_ids[i]
    intact_url += '?content-type=application/json;feature=gene'
    print(intact_url)
    temp = pd.read_json(intact_url)
    for j in range(len(temp.external_name)):
        if pd.isnull(temp.external_name[j]):
            pass
        else:
            genes_and_overlap_node_name.append(temp.external_name[j])
            genes_and_overlap_node_type.append('associated_and_overlapping_genes')
            genes_and_overlap_superphenotype.append(df2_transposed[i][3].split('; ')[-1])
            genes_and_overlap_mim_number.append(df_geneMap2_transposed[i]['entry.mimNumber'])

genes_and_overlap_df = pd.DataFrame(index=range(len(genes_and_overlap_node_name)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

for i in range(len(genes_and_overlap_df)):
    genes_and_overlap_df['Node_name'][i] = genes_and_overlap_node_name[i]
    genes_and_overlap_df['Node_type'][i] = genes_and_overlap_node_type[i]
    genes_and_overlap_df['Superphenotype'][i] = genes_and_overlap_superphenotype[i]
    genes_and_overlap_df['MIM_number'][i] = genes_and_overlap_mim_number[i]

for i in range(len(mimdf)):
    ensembl_df['Node_name'][i] = df_geneMap2_transposed[i]['entry.phenotypeMapList']['phenotypeMap.ensemblIDs'].split(',')[0]
    ensembl_df['Node_type'][i] = 'phenotypeMap.ensemblIDs'
    ensembl_df['Superphenotype'][i] = df2_transposed[i][3].split('; ')[-1]
    ensembl_df['MIM_number'][i] = df_geneMap2_transposed[i]['entry.mimNumber']

for i in range(len(genes_and_overlap_df)):
    genes_and_overlap_df['Superphenotype'][i] = genes_and_overlap_superphenotype[i]
    genes_and_overlap_df['Node_name'][i] = genes_and_overlap_node_name[i]
    genes_and_overlap_df['Node_type'][i] = genes_and_overlap_node_type[i]
    genes_and_overlap_df['MIM_number'][i] = genes_and_overlap_mim_number[i]

gpn = pd.concat([gpn, genes_and_overlap_df], ignore_index=True)


# -------------------------------------------------------------------------------
# First neighbors |
# ----------------+                         (but will remain here commented out)

#intact_url2 = "https://www.ebi.ac.uk/intact/ws/interactor/findInteractor/col25a"
import json
import requests
import ast
unique_list = []
superphenotype_list = []
mim_list = []
for i in range(len(genes_and_overlap_df)):

    intact_url2 = 'https://www.ebi.ac.uk/intact/ws/interactor/findInteractor/'
    intact_url2 += str(genes_and_overlap_df['Node_name'][i])

    df_new = json.loads(requests.get(intact_url2).text)
    array = {}
    for j in range(len(df_new['content'])):
        array[j] = df_new['content'][j]['interactorName']
    my_list = ",".join(array.values()).split(',')
    # traverse for all elements
    for j in my_list:
        superphenotype_list.append(genes_and_overlap_df['Superphenotype'][i])
        mim_list.append(genes_and_overlap_df['MIM_number'][i])
        # check if exists in unique_list or not
        if j not in unique_list:
            unique_list.append(j)
# -------------------------------------------------------------------------------
neighbors_df = pd.DataFrame(index=range(len(unique_list)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

for i in range(len(unique_list)):
    neighbors_df['Node_name'][i] = unique_list[i]
    neighbors_df['Node_type'][i] = 'Intact_first_neighbors'
    neighbors_df['Superphenotype'][i] = superphenotype_list[i]
    neighbors_df['MIM_number'][i] = mim_list[i]

gpn = pd.concat([gpn,neighbors_df],ignore_index=True)
gpn.drop(columns = 'Node_name_temp', inplace = True)
print(neighbors_df)
"""

