

import argparse
import os
import urllib.parse
import urllib.request



# Parse command line input and options
import pandas as pd

parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Convert your list of gene IDs! ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_LIST_FILENAME.txt>')
parser.add_argument('-s', '--symbol', type=str, help='ENSG, EMBL, etc. (see README_conversion_IDs_list.txt for full list)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_LIST_FILENAME>')
args = parser.parse_args()

input = args.input
output = args.output
symbol = args.symbol

logo = """
 +-------------------------------------------------------++     O=o
 |/ o --  __            __    __     ~ ` `       ---<  o ||      O   C l e m s o n
 |/      / _  _ _  _  ||  \  /   _  _    _ _|_ _ ,_     _||     o=O  
 |/_     \__)(-| )(-  ||__/  \__(_)| )\/(-| |_(-|        ||    0===0     U n i v e r s i t y   
 |/__              _   __                            _ __||     O=o       
 |/ o   _____    '        --`       ₲Ɇ₦Ø₱ⱧɆ₦Ø  v5.4    o ||      O          |\_/|
 |/______________________________________________________||      o=O      =( o O )=
 |┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴┬┴|    0===0       /\ " /\\
                                                                O=o       | |\_/| |
                                                                          | |\_/| |
                                                                          \_>---<_/
                                                                          (___|___)
"""
print(logo)


# opening the file in read mode
my_file = open(input, "r")

# reading the file
data = my_file.read()

# replacing end splitting the text
# when newline ('\n') is seen.
input_list = data.split("\n")

my_file.close()

query_string = ''
for i in range(len(input_list)):
    query_string += '\"'
    query_string += input_list[i]
    query_string += '\",'
    query_string += ' '

import requests
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
    json={
        'organism':'hsapiens',
        'target':symbol,
        'query':input_list,
    }
    )
result = r.json()['result']

df = pd.DataFrame(result)

incoming_id_list = []
incoming_id_list2 = []
for i in range(len(df)):
    if df['incoming'][i] == 'None':
        pass
    else:
        incoming_id_list += [df['incoming'][i]]



for i in incoming_id_list:
    if i == '':
        pass
    else:
        incoming_id_list2 += [i]

print('Your input list:')
print()
k=0
print(len(incoming_id_list))


for i in range(len(incoming_id_list2)):
    if i < 6:
        print(incoming_id_list2[i])
        k += 1
    else:
        break


if k > 5 & k < len(incoming_id_list2):
    print('And so on...')
print()





bad_id_list = []

textfile = open(output, "w")

for i in range(len(df)):
    if df['converted'][i] == 'None':
        bad_id_list += [df['incoming'][i]]
    else:
        textfile.write(df['converted'][i] + "\n")

textfile.close()

bad_id_output_filename = output.split('.')[0]+'_FAILED_IDs.txt'

textfile = open(bad_id_output_filename, "w")

for i in bad_id_list:
    if i == '':
        pass
    else:
        textfile.write(i + "\n")
textfile.close()

print(' * Converted your list of',len(input_list),'gene IDs to',symbol+'.')
print()
print(' *',len(df)-len(bad_id_list),'of',len(df),'IDs were successfully converted and saved to \"'+output+'\".')
print()
print(' * A list of the',len(bad_id_list),'failed IDs was saved to \"'+bad_id_output_filename+'\".')
print()