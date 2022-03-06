import numpy as np
import pandas as pd
import argparse
import json
import requests
import ast
import matplotlib.pyplot as plt
from io import BytesIO

logo = """
      _______           __   ____      __                       __                 
     / ____(_)___  ____/ /  /  _/___  / /____  _________ ______/ /_____  __________
    / /_  / / __ \/ __  /   / // __ \/ __/ _ \/ ___/ __ `/ ___/ __/ __ \/ ___/ ___/
   / __/ / / / / / /_/ /  _/ // / / / /_/  __/ /  / /_/ / /__/ /_/ /_/ / /  (__  ) 
  /_/   /_/_/ /_/\__,_/  /___/_/ /_/\__/\___/_/   \__,_/\___/\__/\____/_/  /____/  

                                                        [ G e n o P h e n o ]  3.2                         
"""
print(logo)


# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Find first interactors for the genes in your genotype/phenotype table! * ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_FILENAME.csv>  (phenotype table with ENSEMBL IDs)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME.csv>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output

gpn = pd.read_csv(input)

# Drop the unwanted default index that arises after using the pd.concat method :(
gpn.drop(columns='Unnamed: 0',inplace=True)
print()
print('--------------------------------------------------------------------')
print('	ʕっ•ᴥ•ʔっ        Your input table...   ')
print()
print()
print(gpn)
print()
print()
print('Extracting approved gene symbols from input table:')
print()
gene_ids_list = []
gene_ids_list_unique = []

for i in range(len(gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols'].reset_index())):
    gene_ids_list += [gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols'].reset_index()['Node_name'][i]]

for i in gene_ids_list:
    if i in gene_ids_list_unique:
        pass
    else:
        gene_ids_list_unique += [i]

ensembl_ids_list = []
ensembl_ids_list_unique = []

for i in range(len(gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs'].reset_index())):
    ensembl_ids_list += [gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs'].reset_index()['Node_name'][i]]

for i in ensembl_ids_list:
    if i in ensembl_ids_list_unique:
        pass
    else:
        ensembl_ids_list_unique += [i]

query_string = ''
print(ensembl_ids_list_unique)
for i in ensembl_ids_list_unique:
    query_string += i
    query_string += ' '

import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ENSEMBL_ID',
'to': 'ACC',
'format': 'tab',
'query': query_string
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()

conversion = response.decode('utf-8')
conversion2 = conversion.split('\n')

protein_map_temp = []
ensembl_map_temp = []
protein_map = []
ensembl_map = []

for i in conversion2:
    ensembl_map_temp += [i.split('\t')[0]]
    if '\t' in i:
        protein_map_temp += [i.split('\t')[1]]

for i in range(len(protein_map_temp)):
    if len(protein_map_temp[i]) == 6:
        protein_map += [protein_map_temp[i]]

intact_url = ''
dictionary = {}
df = pd.DataFrame()

for i in range(len(ensembl_ids_list_unique)):
    intact_url = 'https://www.ebi.ac.uk/intact/ws/interaction/list?intraSpeciesFilter=true&draw=50&maxMIScore=1&minMIScore=0&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000&query='
    intact_url += protein_map[i]
    r = requests.post(intact_url)
    response = r.json()
    dictionary.update(response)
    tempdf = pd.DataFrame.from_dict(dictionary['data'])
    if i == 0:
        df = pd.DataFrame.from_dict(dictionary['data'])
    else:
        df = pd.concat([df,tempdf],axis=0,ignore_index=True)

intact_url = ''
dictionary = {}
df2 = pd.DataFrame()
for i in range(len(ensembl_ids_list_unique)):
    intact_url = 'https://www.ebi.ac.uk/intact/ws/interaction/list?intraSpeciesFilter=true&draw=50&maxMIScore=1&minMIScore=0&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000&query='
    intact_url += ensembl_ids_list_unique[i]
    r = requests.post(intact_url)
    response = r.json()
    dictionary.update(response)
    tempdf = pd.DataFrame.from_dict(dictionary['data'])
    if i == 0:
        df2 = pd.DataFrame.from_dict(dictionary['data'])
    else:
        df2 = pd.concat([df2, tempdf], axis=0, ignore_index=True)

    print(i, 'OMIM Gene:         ', protein_map[i], '             ENSMBL ID:         ', ensembl_ids_list_unique[i],           '              Total Rows:          ', len(tempdf), 'df rows:  ', len(df2) + len(df))

df = pd.concat([df, df2], axis=0, ignore_index=True)


"""
# Potentially useful code for the future, leftover from migration from table.py
    for j in range(len(temp_df)):
        df.loc[k] = [temp_df['moleculeA'][j], temp_df['moleculeB'][j]]
        k += 1
"""




"""
# -------------------------------------------------------------------------------------------------
# Intact Material |
#-----------------+

ensembl_df = pd.DataFrame(columns=['Superphenotype', 'Node_name', 'Node_type', 'MIM_number','Parenthetical','Node_name_temp'])
gpn = pd.read_csv(args.input)

# Drop the unwanted default index that arises after using the pd.concat method :(
gpn.drop(columns='Unnamed: 0',inplace=True)
print()
print('--------------------------------------------------------------------')
print('	ʕっ•ᴥ•ʔっ        Your input table...   ')
print()
print()
print(gpn)
print()

gene_ids = []
for i in gpn[gpn['Node_type'] == 'phenotypeMap.approvedGeneSymbols']['Node_name']:
    gene_ids += [i]
    print(i)
print()
print('--------------------------------------------------')
print(len(gene_ids),' Gene IDs imported from table')
print('--------------------------------------------------')


print()
print('Extracting ENSMBL IDs...')
print()

ensembl_ids = []
for i in gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs']['Node_name']:
    ensembl_ids += [i]
    print(i)
print()
print('--------------------------------------------------')
print(len(ensembl_ids),' ENSMBL IDs imported from table')
print('--------------------------------------------------')

print()
print('Finding overlapping genes for given ENSMBL IDs:')
print()
intact_url = ''
genes_and_overlap = []
for i in range(len(ensembl_ids)):
    intact_url = 'https://rest.ensembl.org/overlap/id/'
    intact_url += ensembl_ids[i]
    intact_url += '?content-type=application/json;feature=gene'
    temp = pd.read_json(intact_url)
    for j in range(len(temp.external_name)):
        if pd.isnull(temp.external_name[j]):
            pass
        else:
            genes_and_overlap.append(temp.external_name[j])
            print('Overlap found:              ',ensembl_ids[i])

print()
print('---------------------+')
print('Overlapping features |')
print('---------------------+-----------------------------------------------------------------------------------')
print(genes_and_overlap)
print('---------------------------------------------------------------------------------------------------------')
print()

# For testing: (remove later)
#for i in genes_and_overlap:
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
# First neighbors | ** ADD CODE HERE **
# ----------------+

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
print(df)
# Save the gpn as a csv using the same filename, but with extension '.csv'
df.to_csv(output)

