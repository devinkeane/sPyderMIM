# Import libraries
import argparse
import os

def save_to_file(text,output):
    with open(output, mode='wt', encoding='utf-8') as myfile:
        myfile.write('\n'.join(text))
        myfile.write('\n')
# ------------------------------------------------------------------------------------------------------
# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Split one large MIM list into lists of 20 or less MIMs to use in table.py (for functionality, due to API limits)  * ")
parser.add_argument('-i', '--input', type=str, help='<INPUT_FILENAME.txt>   (One MIM list larger than 20 MIMs to be divided)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME>    (without \".txt\":  multiple numbered .txt files will be produced automatically.)')
args = parser.parse_args()

input = args.input
output = args.output


opening_screen = """

   _       ,/'
  (_).  ,/'
   _  ::       ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  
  (_)'  `\.
           `\.              [  G e n o P h e n o  ]
                                        List Splitter
"""
print()
print(opening_screen)
print()


my_file = open(input, "r")
lst = my_file.read().splitlines()

chunked_list = []
chunk_size = 20

for i in range(0, len(lst), chunk_size):
    chunked_list.append(lst[i:i+chunk_size])

dictionary = {}
print('Saving divided lists to directory:    '+str('./{}_MIM_directory'.format(output))+ str('/'))
print()

newpath = './{}_MIM_directory'.format(output)
if not os.path.exists(newpath):
    os.makedirs(newpath)

for i in range(len(chunked_list)):
    dictionary['./{}_MIM_directory/'.format(output)+output+"{0}.txt".format(i)] = chunked_list[i]
    print(str(output) + str(i) + '.txt')
print()
for i in dictionary:
    save_to_file(dictionary[i], i)