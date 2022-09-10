# last rev 2022-09-07

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


# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Find first interactors for the genes in your genotype/phenotype table! * ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_FILENAME.csv>  (phenotype table with ENSEMBL IDs)')
parser.add_argument('-m', '--mode', type=str, help='\"omim\" (genes or clinical data output from table.py) or \"list\" (simple .txt list of ENSG IDs)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME.csv>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
mode = args.mode

# -------------------------------------------------------------------------------------------
# Loading bar animation function
done = False
def animate():


    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            sys.stdout.write('')
            sys.stdout.flush()
            break
        sys.stdout.write('\r(⌐ ͡■ ͜ʖ ͡■) Working on it... ' + c)
        sys.stdout.flush()
        sys.stdout.write('\r')
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.flush()
    sys.stdout.write('\r                           ')
    sys.stdout.write('\r( ͡° ͜ʖ ͡°)ﾉ⌐■-■  We did it! ✔')
    sys.stdout.flush()

searching_wait_animation = threading.Thread(target=animate)

# -------------------------------------------------------------------------------------------

# create and print logo
logo = """
      _______           __   ____      __                       __                 
     / ____(_)___  ____/ /  /  _/___  / /____  _________ ______/ /_____  __________
    / /_  / / __ \/ __  /   / // __ \/ __/ _ \/ ___/ __ `/ ___/ __/ __ \/ ___/ ___/
   / __/ / / / / / /_/ /  _/ // / / / /_/  __/ /  / /_/ / /__/ /_/ /_/ / /  (__  ) 
  /_/   /_/_/ /_/\__,_/  /___/_/ /_/\__/\___/_/   \__,_/\___/\__/\____/_/  /____/  

                                                       [ G e n o P h e n o ]  v7.2                         
"""
print(logo)
# -------------------------------------------------------------------------------------------

# instantiate a dataframe
gpn = pd.DataFrame()


# import genes, either from an OMIM table for from a simple .txt list

if mode == 'omim':
    gpn = pd.read_csv(input)

    # Drop the unwanted default index that arises after using the pd.concat method :(
    gpn.drop(columns='Unnamed: 0',inplace=True)

if mode == 'list':

    bashCommand = 'python3 ./convert_ids.py -i {} -s ENSG -o converted_ENSG.txt >/dev/null 2>&1'.format(input)
    os.system(bashCommand)


    bashCommand = 'python3 ./convert_ids.py -i {} -s HGNC -o converted_HGNC.txt >/dev/null 2>&1'.format(input)
    os.system(bashCommand)

    # opening the file in read mode
    my_file = open('converted_ENSG.txt', "r")

    # reading the file
    data = my_file.read()

    # replacing end splitting the text when newline ('\n') is seen.
    ENSG_input_list = data.split("\n")

    my_file.close()

    my_file = open('converted_HGNC.txt', "r")

    # reading the file
    data = my_file.read()

    # replacing end splitting the text
    # when newline ('\n') is seen.
    HGNC_input_list = data.split("\n")

    my_file.close()

    gpn = pd.DataFrame(index=range(len(ENSG_input_list)), columns=['Node_name', 'Node_type', 'Gene_MIM_number','Phenotype_MIM_number'])
    gpn_ENSG = pd.DataFrame(index=range(len(ENSG_input_list)),columns=['Node_name','Node_type', 'Gene_MIM_number','Phenotype_MIM_number'])
    gpn_HGNC = pd.DataFrame(index=range(len(ENSG_input_list)), columns=['Node_name', 'Node_type', 'Gene_MIM_number','Phenotype_MIM_number'])

    for i in range(len(ENSG_input_list)):
        gpn_ENSG['Node_name'][i] = ENSG_input_list[i]
        gpn_ENSG['Node_type'][i] = 'phenotypeMap.ensemblIDs'

    for i in range(len(HGNC_input_list)):
        gpn_HGNC['Node_name'][i] = HGNC_input_list[i]
        gpn_HGNC['Node_type'][i] = 'phenotypeMap.approvedGeneSymbols'

    gpn = pd.concat([gpn_ENSG,gpn_HGNC])

    my_file = open('converted_ENSG_FAILED_IDs.txt', "r")

    # reading the file
    data = my_file.read()

    # replacing end splitting the text
    # when newline ('\n') is seen.
    bad_id_list = data.split("\n")

    my_file.close()

    bashCommand = 'mv ./converted_ENSG_FAILED_IDs.txt ./{}_FAILED_IDs.txt >/dev/null 2>&1'.format(output.split('.')[0])
    os.system(bashCommand)

    bashCommand = 'rm ./converted_HGNC_FAILED_IDs.txt >/dev/null 2>&1'
    os.system(bashCommand)

    bashCommand = 'rm ./converted_HGNC.txt >/dev/null 2>&1'
    os.system(bashCommand)

    bashCommand = 'rm ./converted_ENSG.txt >/dev/null 2>&1'
    os.system(bashCommand)

print()
print('--------------------------------------------------------------------')
print()
print('	ʕっ•ᴥ•ʔっ        Your input:   ')
print()
print()

if mode == 'omim':
    time.sleep(2)
    print(gpn)

if mode == 'list':
    my_file = open(input, "r")

    # reading the file
    data = my_file.read()
    # replacing end splitting the text
    # when newline ('\n') is seen.
    input_list = data.split("\n")

    if len(input_list) < 11:
        for i in input_list:
            print(i)
    else:
        k = 0
        for i in input_list:
            if k < 11:
                print(i)
                k += 1
            else:
                break
    print('And so on...    *:･ﾟ✧')


    my_file.close()
print()
if mode == 'omim':
    print('Extracting Gene MIM IDs from input:')
if mode == 'list':
    print('Converting input to ENSEMBL IDs:')
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
ensembl_ids_list_unique = list(filter(None, ensembl_ids_list_unique))


gene_MIM_list = []
gene_MIM_list_unique = []

if mode == 'omim':
    for i in range(len(gpn)):
        gene_MIM_list += [gpn['Gene_MIM_number'][i]]

    for i in gene_MIM_list:
        if i in gene_MIM_list_unique:
            pass
        else:
            gene_MIM_list_unique += [i]


query_string = ''

if mode == 'list':
    for i in range(len(ensembl_ids_list_unique)):
        if i < len(ensembl_ids_list_unique) - 1:
            query_string += '\''+ensembl_ids_list_unique[i]+'\','
        else:
            query_string += '\'' + ensembl_ids_list_unique[i] + '\''

    r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism': 'hsapiens',
        'target': 'MIM_GENE_ACC',
        'query': ensembl_ids_list_unique
    }
    )
    result = r.json()['result']

    MIM_conversion_df = pd.DataFrame(result)

    for i in MIM_conversion_df['converted']:
        gene_MIM_list_unique += [i]

query_string = ''

sleep_float = 0.1
if len(ensembl_ids_list_unique) < 1000:
    sleep_float = 0.0
if len(ensembl_ids_list_unique) < 200:
    sleep_float = 0.025
if len(ensembl_ids_list_unique) < 50:
    sleep_float = 0.1
if len(ensembl_ids_list_unique) < 20:
    sleep_float = 0.15
if len(ensembl_ids_list_unique) < 10:
    sleep_float = 0.2

for i in gene_MIM_list_unique:
    sys.stdout.write('\r'+str(i))
    time.sleep(sleep_float)
    sys.stdout.flush()

sys.stdout.write('\r ')
time.sleep(0.2)
sys.stdout.flush()
sys.stdout.write('\r'+'Gene MIMs extracted     ✔       '+'\n')
sys.stdout.flush()
print()

sys.stdout.write('Converting MIM HGNC Gene IDs to UNIPROT IDs       ')
sys.stdout.flush()


query_string = '\''
for i in range(len(gene_MIM_list_unique)):
    query_string += str(gene_MIM_list_unique[i])
    if i < len(gene_MIM_list_unique)-1:
        query_string += ', '

params = {
    'from': (None, 'MIM'),
    'to': (None, 'UniProtKB'),
    'ids': (None, query_string),
}

response = requests.Response()
while response.status_code != 200:
    response = requests.post('https://rest.uniprot.org/idmapping/run', params=params)
job_ID = response.json()['jobId']


response2 = requests.get('https://rest.uniprot.org/idmapping/status/'+job_ID)
while 'results' not in response2.json():
    response2 = requests.get('https://rest.uniprot.org/idmapping/status/'+job_ID)
response2 = requests.get(f'https://rest.uniprot.org/idmapping/results/'+job_ID+'/?size=500')

sys.stdout.write('\r'+'Conversion to UNIPROT complete  ✔                         '+'\n')
time.sleep(2)
sys.stdout.flush()



intact_url = ''
dictionary = {}
df2 = pd.DataFrame()
tempdf = pd.DataFrame()


# check for failed IDs and drop them from their respective mapped lists
if 'failedIds' in response2.json():
    for i in range(len(response2.json()['failedIds'])):
        j = 0
        for x in gene_MIM_list_unique[:]:
            if str(x) == response2.json()['failedIds'][i]:
                gene_MIM_list_unique.pop(j)
                ensembl_ids_list_unique.pop(j)
                gene_ids_list_unique.pop(j)
                j -= 1
            j += 1
        print()
        print('⚠ WARNING ⚠    Failed to convert Gene MIM ' + str(response2.json()['failedIds'][i]) + ' to UNIPROT ID')
        print()


# Query each ID using the UNIPROT ID conversions in the response.

# Print to output each respective ENSEMBL ID and HGNC symbol for each iterations,
# using the same symbols when the next iteration contains another protein transcribed
# by the same gene.
print()
print('                                         ┌(◉ ͜ʖ◉)つ')
print('         +-----------------------------------------------------------------------+')
print('         |   Searching for proteins that interact with each OMIM gene product:   |')
print('         +-----------------------------------------------------------------------+')
print()

# Begin waiting animation
searching_wait_animation.start()

i = 0
for j in range(len(response2.json()['results'])):

    # ---------------------------------------------------------------------
    intact_url = 'https://www.ebi.ac.uk/intact/ws/interaction/list?draw=50&interactorSpeciesFilter=Homo%20sapiens&interactorTypesFilter=protein&intraSpeciesFilter=true&maxMIScore=1&minMIScore=0&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000&query='
    intact_url += response2.json()['results'][j]['to']
    r = requests.post(intact_url)
    while r.status_code == 500:
        r = requests.post(intact_url)
    response = r.json()
    dictionary.update(response)
    tempdf = pd.DataFrame.from_dict(dictionary['data'])


    r = requests.post(intact_url)
    while r.status_code == 500:
        r = requests.post(intact_url)
    response = r.json()
    dictionary.update(response)
    tempdf = pd.DataFrame.from_dict(dictionary['data'])

    if j == 0:
        df2 = pd.DataFrame.from_dict(dictionary['data'])
    else:
        df2 = pd.concat([df2, tempdf], axis=0, ignore_index=True)
    sys.stdout.flush()

    if j > 0:
        if response2.json()['results'][j]['from'] == response2.json()['results'][j-1]['from']:
            i = i - 1

    print(j+1, '| Approved Gene ID:', gene_ids_list_unique[i], '| Gene Product (UNIPROT ID):', response2.json()['results'][j]['to'], '| Interactions:', len(tempdf), '| Total Interactions:', len(df2))
    sys.stdout.flush()
    i = i+1


print()

done = True

time.sleep(2)

#df = pd.concat([df, df2], axis=0, ignore_index=True)

sys.stdout.flush()
print()
print()
print('Saving output table:')
print()
print(df2)
print()
print()
print('Process complete.')
print()

if mode == 'list':
    if len(bad_id_list) > 1:
        print('-------------------------------------------------------------------------------------------------------------')
        print()
        print()
        print('  *  '+str(len(bad_id_list)-1)+' of your original input IDs failed to undergo conversion to ENSEMBL IDs, which interactors.py requires.')
        print()
        print('           *  These IDs have been saved to \"{}_FAILED_IDs.txt\".'.format(output.split('.')[0]))
print()
print('-------------------------------------------------------------------------------------------------------------')
print()
print('Your interactors list was saved as \"'+output+'\" with ',len(df2.columns),' columns and ',len(df2),' rows.')
print()
print('                                                    	⊂(◉‿◉)つ  ♡  ')

print()
# Save the gpn as a csv
df2.to_csv(output)

# END OF PROGRAM





# - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . -
#                                                                                              ;
# Potentially useful code for the future is below, leftover from various phases of devolopment ;
# This code will likely be in favor of the evolving Master's Thesis plan                       ;
#           -DK  2022-09-10                                                                    ;
#                                                                                              ;
# - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . - , - . - . - . -

# ---------------------------------------------------------------------

"""
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism': 'hsapiens',
        'target': 'UNIPROT_GN_ACC',
        'query': ensembl_ids_list_unique[i],
    }
)

result = r.json()['result']


uniprot_conversion_df = pd.DataFrame(result)
"""


# The following portion of code is likely to be deprecated in the near future
# G:Profiler now performs this task instead
"""
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
    if '\t' in i:
        ensembl_map_temp += [i.split('\t')[0]]
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
"""


"""


  
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
print('Extracting ENSEMBL IDs...')
print()

ensembl_ids = []
for i in gpn[gpn['Node_type'] == 'phenotypeMap.ensemblIDs']['Node_name']:
    ensembl_ids += [i]
    print(i)
print()
print('--------------------------------------------------')
print(len(ensembl_ids),' ENSEMBL IDs imported from table')
print('--------------------------------------------------')

print()
print('Finding overlapping genes for given ENSEMBL IDs:')
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