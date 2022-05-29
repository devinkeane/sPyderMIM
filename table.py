#!/usr/bin/env python
# coding: utf-8

#     ___|                       _ \   |
#    |       _ \  __ \    _ \   |   |  __ \    _ \  __ \    _ \   |
#    |   |   __/  |   |  (   |  ___/   | | |   __/  |   |  (   |  |
#   \____| \___| _|  _| \___/  _|     _| |_| \___| _|  _| \___/   |
#   ______________________________________________________________|
#                                      (ツ)_/¯  - * Version 5.2 * -
#  [ O m i m   T a b l e   M a k e r ]
#
# Last rev: 2022-03-12
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
                                      (ツ)_/¯  - * Version 5.5 * -
  [ O M I M   T a b l e   M a k e r ]
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
    print()
    print('Please limit your MIM list to 20 numbers per \".txt\" file.')
    print('You can then use \"concat.py\" to combine each of these outputs')
    print('into one table if you need to query more than 20 MIMs.')
    print()
    exit(-1)

# ------------------------------------------------------------------------------------------------------
# Inserting the MIM numbers into the API request URL and retrieving data |
# -----------------------------------------------------------------------+

# Begin building the OMIM API request URL
url = 'https://api.omim.org/api/entry?mimNumber='

# Append the list of MIM numbers (now in "mimdf") onto URL
for line in range(len(mimdf)):
    if type(mimdf[0][line]) == str:
        url += mimdf[0][line]
        url += ','
    else:
        url += mimdf[0][line].astype(str)
        url += ','

# Append clinical synopsis API request, JSON format option onto URL
url += '&include=clinicalSynopsis&format=json&apiKey='

# Append API key from parsed command line arguments onto URL
url += args.apikey

# Create a second URL for the gene map API, simply by replacing 'clinicalSynopsis' with 'geneMap'
# at the the '&include' statement in the API URL
url_geneMap = url.replace('clinicalSynopsis','geneMap')

# ------------------------------------------------------------------------------------------------------
# Making the API requests and storing the data in a dataframe  |
# -------------------------------------------------------------+


# Load data into new data frame using pandas .read_json() module
df = pd.read_json(url)

# flatten data (data is entirely nested at one column 'omim' and one row @ index 0)
df2 = pd.json_normalize(df['omim'][0])
# - - - - - - - - - -
# load geneMap data using Python JSON module
df_geneMap = pd.read_json(url_geneMap)

# This dataframe may be renamed in the future.  It ends up being the final dataframe that
# exported, but it is no longer transposed as a result of changes to the code over time.
df_geneMap2_transposed = pd.DataFrame()

# flatten data
for i in range(len(pd.json_normalize(df_geneMap['omim'][0]))):
    df_geneMap2 = pd.json_normalize(df_geneMap['omim'][0][i])
    temp = pd.DataFrame.transpose(df_geneMap2)
    df_geneMap2_transposed = pd.concat([df_geneMap2_transposed, temp], axis=1,ignore_index=True)

# drop molecular basis
# (not so data friendly gene id column, we will use the nested geneMap list from the
# API request instead)
if 'entry.clinicalSynopsis.molecularBasis' in df2.columns:
    df2.drop(columns='entry.clinicalSynopsis.molecularBasis',inplace=True)

# Transpose the clinical data --> now each column represents a disease (MIM#)
df2_transposed = pd.DataFrame.transpose(df2)


# The following section may be changed or removed.  It was put in place to keep
# track of failed MIMs to print to output.  However, previous changes I made to
# parsing appear have fixed the program so that it does not ever seem to fail on
# any phenotypic MIM (those which begin with the '#' symbol).  It is worth noting
# that this part of the code appeared to be imperfect anyway and would need to
# be modified if I decide to keep  it.

bad_mim_count = 0

for i in df2_transposed.columns:
    if isinstance(df2_transposed[i]['entry.clinicalSynopsis.inheritance'],float):
        print('No clinical data found...        MIM # :    ',df2_transposed[i]['entry.mimNumber'])
        bad_mim_count += 1

# ----------------------------------------------------------------------------------------------+
# Populating a fresh genotype/phenotype dataframe with data from df2_transposed (OMIM API data) |
# ----------------------------------------------------------------------------------------------+

# Create a new data frame that will be the output and call it "gpn" (genotype/phenotype network)
gpn = pd.DataFrame(columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

empty_row_df = pd.DataFrame([[np.nan] * len(gpn.columns)],columns=gpn.columns)

# Reset indices:  go back to the top row of gpn and row at index 5 on df2_transposed


k = 0
i = 5

# for every row of every column in df2_transposed, starting on row = 5...
for j in range(len(df2_transposed.columns)):
    for row in range(len(df2_transposed)):
        if row > 4:

            # If the current cell is null/float/integer, skip it
            if pd.isnull(df2_transposed[j][row]) | pd.api.types.is_float(df2_transposed[j][row]) | pd.api.types.is_integer(df2_transposed[j][row]):
                 pass

            # If the cell is not null but no '\n' is found, there is only one element in the cell.
            # Add this element to the 'Node_name_temp' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') == 0:
                tempdf = empty_row_df
                tempdf['Node_name_temp'] = df2_transposed[j][row]
                tempdf['Node_type'] = df2_transposed.index[row]
                tempdf['MIM_number'] = df_geneMap2_transposed[j]['entry.mimNumber']

                tempdf['Superphenotype'] = df2_transposed[j][3].split('; ')[-1]

                gpn = pd.concat([gpn, tempdf], ignore_index=True)

            # Otherwise, the amount of elements should be equal to n + 1,
            # so add each of these elements to the 'Node_name_temp' column of gpn at row k.
            elif df2_transposed[j][row].count('\n') > 0:
                for element in range(len(df2_transposed[j][row].split('\n'))):
                    tempdf = empty_row_df
                    tempdf['Node_name_temp'] = df2_transposed[j][row].split('\n')[element].replace(';','')
                    tempdf['Node_type'] = df2_transposed.index[row]
                    tempdf['MIM_number'] = df_geneMap2_transposed[j]['entry.mimNumber']


                    tempdf['Superphenotype'] = df2_transposed[j][3].split('; ')[-1]
                    gpn = pd.concat([gpn, tempdf], ignore_index=True)
                    k += 1


# ---------------------------------------------------------------------------
# Now, we need to clean up the Node_name_temp column data in gpn by cutting the
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

# Create permanent 'Node_name' column from 'Node_name_temp, but without the parentheticals
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

# Instantiate a new dataframe to populate with ENSEMBL ID info   ** MARKED FOR POSSIBLE REMOVAL **
# ensembl_df = pd.DataFrame(index=range(len(mimdf)),columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])


# Parsing ENSEMBL IDs from df_geneMap2_transposed and other data from df2_transposed to write into gpn2
gpn2 = pd.DataFrame(index=gpn,columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])

bad_mim_count2 = 0
print()
k = 0
for i in range(len(df_geneMap2_transposed.columns)):
    # Due to increased accuracy of the program, I intend to depracate this if branch.
    if isinstance(df_geneMap2_transposed[i]['entry.phenotypeMapList'],float):
        print('No ENSEMBL/Gene ID found...      MIM # :    ',df_geneMap2_transposed[i]['entry.mimNumber'])
        bad_mim_count2 += 1
    else:
        for j in range(len(df_geneMap2_transposed[i]['entry.phenotypeMapList'])):
            gpn2['Node_name'][k] = df_geneMap2_transposed[i]['entry.phenotypeMapList'][j]['phenotypeMap']['ensemblIDs'].split(',')[0]
            gpn2['Node_type'][k] = 'phenotypeMap.ensemblIDs'
            gpn2['Superphenotype'][k] = df2_transposed[i][3].split('; ')[-1]
            gpn2['MIM_number'][k] = df_geneMap2_transposed[i]['entry.mimNumber']
            k += 1

# Parsing HUGO gene symbol IDs from df_geneMap2_transposed and other data from df2_transposed to write into gpn2
gpn3 = pd.DataFrame(index=gpn,columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])
k = 0
for i in range(len(df_geneMap2_transposed.columns)):
    if isinstance(df_geneMap2_transposed[i]['entry.phenotypeMapList'],float):
        pass
    else:
        for j in range(len(df_geneMap2_transposed[i]['entry.phenotypeMapList'])):
            gpn3['Node_name'][k] = df_geneMap2_transposed[i]['entry.phenotypeMapList'][j]['phenotypeMap']['approvedGeneSymbols']
            gpn3['Node_type'][k] = 'phenotypeMap.approvedGeneSymbols'
            gpn3['Superphenotype'][k] = df2_transposed[i][3].split('; ')[-1]
            gpn3['MIM_number'][k] = df_geneMap2_transposed[i]['entry.mimNumber']
            k += 1
# Concatenate by appending ENSEMBL IDs (gpn2) and HUGO symbols (gpn3) to original data frame (gpn)
gpn = pd.concat([gpn,gpn2], ignore_index=True)
gpn = pd.concat([gpn,gpn3], ignore_index=True)

# -------------------------------------------------+
# Printing output and saving the dataframe to file |
# -------------------------------------------------+

# Create a list of ENSEMBL IDs from GPN to use for output reporting.
ensembl_ids = (gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs'])['Node_name'].reset_index(drop=True)

gpn.drop(columns='Node_name_temp',inplace=True)
gpn.dropna(how='all',inplace=True)
gpn.reset_index(inplace=True,drop=True)

# Print a preview of gpn to output
print()
print(gpn)
print()

# Save the gpn as a csv using the same filename, but with extension '.csv'
gpn.to_csv(output)

#df2_transposed.to_csv('testing.csv')

# Print output message
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print('  ...Your network table was saved as \"',output,'\"')
print('  with ',len(gpn),' total rows from your list of ',len(mimdf),' MIMs.')
print()
print('  ',len(mimdf)-bad_mim_count,' of ',len(mimdf),' MIMs contained phenotypic data.')
print()
print('  ',len(mimdf)-bad_mim_count2,' of ',len(mimdf),' MIMs contained ENSEMBL/Gene IDs.')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
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
# Intact Material | ** This section is being transferred to interactors.py **
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