 ___   ____   __    ___   _      ____
| |_) | |_   / /\  | | \ | |\/| | |_
|_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__
                               ₲Ɇ₦Ø₱ⱧɆ₦Ø

GenoPheno is a genotype/phenotype network generator that uses the OMIM (Online Mendelian
Inheritance in Man) database in order to construct a relationship graph that links diseases
by nodes that represent genes and phenotypes.

GenoPheno takes in a .txt file that is a list of "MIM numbers," each of which correspond to a
disease listed in the OMIM Database.  In the current version, the output is a .png of a graph
where the source nodes are diseases listed in OMIM.  Additionally, it exports a .csv table which
can be used for further graph theory analysis in other programs.  However, the future goal in
development of this program is to provide a purely command line interface that can allow for
graph theory analysis to be upscaled and applied to larger data sets than could be handled
in GUI-based programs.

Version 2.0 will include the options to generate graphs using genes, phenotypes, or both as
the source nodes, whereas the current version only exports a map that uses both as the source
node.

** Upcoming changes in 2.0: **

    * I'm going to split GenoPheno.py into multiple python scripts that can be executed with a PBS
    or bash script together.  I will include an example of this workflow in the updated README.txt.

    * One program will be used to generate a table, and a separate program will allow for graph
    analysis on the table produced.

        * This will allow me to add an increasing amount of command line arguments and functionality
        that offers more options for graphing the table.

    * After this, I hope to add a script that will concatenate any .csv tables created.

  |   |                           
  |   |   __|   _` |   _` |   _ \ 
  |   | \__ \  (   |  (   |   __/ 
  \___/ ____/ \__,_| \__, | \___| 
                     |___/
# ----------------------------------------------------------------------------

1) Obtain an API key through OMIM:
  --> https://www.omim.org/api

# ----------------------------------------------------------------------------

2) Create list of OMIM reference ids ("MIM" numbers) (20 MAXIMUM!):

(i.e.)
~$ vim input_list.txt
------
615291              [input_list.txt]
614505                      |
120580      <---------------`
...
------
# ----------------------------------------------------------------------------

3) Set up an Anaconda environment with the necessary dependencies:

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib networkx
~$ conda activate GenoPheno

# ----------------------------------------------------------------------------

4) Execute with the following syntax/options:

~$ python GenoPheno -i <input_list.txt> -o <output_file.csv> -a <api_key>

# ----------------------------------------------------------------------------

GenoPheno (c) 2022-01-27 by Devin Keane
