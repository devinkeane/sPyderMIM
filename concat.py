#!/usr/bin/python
#
#           C o n c a t . p y |
#           ------------------+
#                 [ GenoPheno ]
#                                      Last rev: 2022-03-07
# ------------------------------------------------------------------------------------------------------
# Import libraries
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------------------------------------
# Parse command line input and options
parser = argparse.ArgumentParser(description="	ʕっ•ᴥ•ʔっ  * Combine your genotype/phenotype and/or interactors tables into a single one for graph analysis! * ")
parser.add_argument('-i', '--input', nargs='+', type=str, help='<INPUT_FILENAME_1.csv> <INPUT_FILENAME_2.csv> ... (two or more files to concatenate)')
parser.add_argument('-o', '--output', type=str, help='<OUTPUT_FILENAME.csv>')
args = parser.parse_args()

# Assign parsed arguments into local variables
input = args.input
output = args.output
# ------------------------------------------------------------------------------------------------------

# Create an empty dataframe
gpn = pd.DataFrame()

# Iterate through the list of input arguments (file names), read each file into a new dataframe,
# and append each of these to the master (gpn) dataframe
for i in range(len(args.input)):
    temp = pd.read_csv(args.input[i])
    gpn = pd.concat([gpn,temp],ignore_index=True)

# Drop the unwanted default index that arises after using the pd.concat method :(
gpn.drop(columns='Unnamed: 0',inplace=True)

# Print final .csv table to output
gpn.to_csv(output)

# ---------------------------------------------------------------------------
# Print output message |
# ---------------------+

print()
print('           c o n c a t . p y |')
print('           ------------------+')
print('                 [ GenoPheno ]')
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
print('                 ...Your combined table was saved as \"',output,'\"')
print('                 with ',len(gpn),' total rows.')
print()
print('--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+')
print()
