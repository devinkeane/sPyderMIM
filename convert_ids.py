import argparse
import os
import urllib.parse
import urllib.request


# Parse command line input and options
import pandas as pd

parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Convert your list of gene IDs! ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_LIST_FILENAME.txt>')
parser.add_argument('-s', '--symbol', type=str, help='ENSG, EMBL, etc. (see README_conversion_IDs_list.TXT for full list)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_LIST_FILENAME>')
args = parser.parse_args()

input = args.input
output = args.output
symbol = args.symbol

logo = """
 +------------------------------------------------------------++     O=o
 |/ o  __            ___ _     __                           o ||      O
 |/   /__ _ __  _     | | \   /   _ __     _  ___|_ _  __   --||     o=O
 |/_  \_|(/_| |(/_   _|_|_/   \__(_)| |\_/(/_ |  |_(/_ |      ||    0===0
 |/__                                                     _ __||     O=o
 |/ o   _____                             ₲Ɇ₦Ø₱ⱧɆ₦Ø  5.2    o ||      O
 |/////////////////////////////////////////////////////////////|     o=O
                                                                   0===0
                                                                    O=o      
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

bad_id_list = []

textfile = open(output, "w")

for i in range(len(df['converted'])):
    if df['converted'][i] == 'None':
        bad_id_list += [df['incoming'][i]]
    else:
        textfile.write(df['converted'][i] + "\n")

textfile.close()

bad_id_output = output.split('.')[0]+'_FAILED_IDs.txt'

textfile = open(bad_id_output, "w")

for i in bad_id_list:
    textfile.write(i + "\n")
textfile.close()

print(' * Converted your list of',len(input_list),'gene IDs to',symbol+'.')
print()
print(' *',len(df)-len(bad_id_list)-1,'of',len(df)-1,'IDs were successfully converted and saved to \"'+output+'\".')
print()
print(' * A list of the',len(bad_id_list),'failed IDs was saved to \"'+bad_id_output+'\".')
print()