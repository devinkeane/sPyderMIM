  ___                       _                             _               _
  | _ \    ___    __ _    __| |   _ __     ___            | |_    __ __   | |_   
  |   /   / -_)  / _` |  / _` |  | '  \   / -_)     _     |  _|   \ \ /   |  _|  
  |_|_\   \___|  \__,_|  \__,_|  |_|_|_|  \___|   _(_)_   _\__|   /_\_\   _\__|  
_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""| 
"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'  
                                                     ₲Ɇ₦Ø₱ⱧɆ₦Ø

GenoPheno is a genotype/phenotype network generator that uses the OMIM database
in order to construct a relationship graph that links diseases by nodes that represent
genes and phenotypes.

The program appears to function properly so far, but this is the beta version first release.
It takes in a .txt file that is a just a list of "mim numbers," each of which correspond to a
disease listed in the OMIM (Online Mendelian Inheritance in Men) Database.  In this version,
the output is a .png of a graph where the source nodes are diseases listed in OMIM.  Additionaly, 
it exports a .csv table which can be used for further graph theory analysis in other programs.

The next version will include the ability to generate graphs using genes, phenotypes, or both as
the source nodes.
              
  |   |                           
  |   |   __|   _` |   _` |   _ \ 
  |   | \__ \  (   |  (   |   __/ 
  \___/ ____/ \__,_| \__, | \___| 
                     |___/

1) Obtain an API key through OMIM:
  --> https://www.omim.org/api

2) Create list of OMIM reference ids ("MIM" numbers):

(i.e.)
~$ vim input_list.txt
------
615291              [input_list.txt]
614505                      |
120580      <---------------`
...
------

3) Set up an Anaconda environment with the necessary dependencies:

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib networkx
~$ conda activate GenoPheno

4) Execute with the following syntax/options:

~$ python GenoPheno -i <input_list.txt> -o <output_file.csv> -a <api_key>

GenoPheno (c) 2022-01-27 by Devin Keane
