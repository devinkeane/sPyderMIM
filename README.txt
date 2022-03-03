 ___   ____   __    ___   _      ____
| |_) | |_   / /\  | | \ | |\/| | |_
|_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__
                               ₲Ɇ₦Ø₱ⱧɆ₦Ø

GenoPheno is a workflow suite and genotype/phenotype network generator that uses the OMIM
(Online Mendelian Inheritance in Man) database in order to construct a relationship graph that
links diseases by genes and phenotypes.  The workflow also utilizes the IntAct API, which
allows the user to find the protein products for the genes associated with each MIM number
and the proteins known to interact with each of this.  Large tables and network graphs can
be constructed in order to investigate dozens of diseases

table.py takes in a .txt file that is a list of "Phenotype MIM numbers," each of which correspond
to a disease subtype listed in the OMIM Database.  A disease subtype in the OMIM database is
defined by a gene that is known to give rise to a set of phenotypes, which are classified
together under the clinical data for that subtype in the OMIM database.  The only required
input for this workflow is a .txt file(s) that is a simple list of phenotype MIM numbers that
the user has collected for diseases they wish to investigate, each separated by a new line.

interactors.py uses the output from table.py to find the protein products of each gene and
their protein interactors.

graph.py is able to take in the .csv output from either table.py or interactors.py and
generate a network graph of phenototypes or protein interactors, respectively.

Output from these tables can be used as edge list input for other programs, such as Cytoscape!
However, the goal of this program is to provide a purely command line interface that can allow
for graph theory analysis to be upscaled and applied to larger data sets than could be handled
in GUI-based programs.

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
** IMPORTANT NOTES FOR VERSION 3.1 ! **   <-- Read first to ensure functionality
---------------------------------------
    ->  GenoPheno.py, the original version of the program, has been replaced by
    table.py, but it will still remain on GitHub for the time being.  table.py
    may later be renamed to GenoPheno.py after the original version is removed.

    ->  GenoPheno currently only takes in lists of "Phenotype MIM numbers" as the
    starting input of the workflow.  DO NOT use "Gene/Locus MIM numbers" as these
    will cause the program to crash.

    ->  When creating a .txt list mim numbers, ONLY USE MIM numbers with
    the "#" prefix.  The current version of this program does not process MIM
    numbers with the "+", "%", or "^" prefixes.

    ->  DO NOT use mim numbers for listings where "susceptibility to" is
    is included in the subtype title.  The current version of this program
    does not support MIM numbers for this or individual genes.

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
  |   |
  |   |   __|   _` |   _` |   _ \
  |   | \__ \  (   |  (   |   __/
  \___/ ____/ \__,_| \__, | \___|
                     |___/
--------------------------------------------------------------------------------------------

1) Obtain an API key through OMIM:
  --> https://www.omim.org/api

--------------------------------------------------------------------------------------------

2) Create list of OMIM reference ids ("MIM" numbers) (20 MAXIMUM!):

(i.e.)
~$ vim input_list.txt
------
615291              [input_list.txt]
614505                      |
120580      <---------------`
...
------
--------------------------------------------------------------------------------------------

3) Set up an Anaconda environment with the necessary dependencies:

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib networkx
~$ conda activate GenoPheno

--------------------------------------------------------------------------------------------

4) Execute with the following syntax/options:
--------------------------------------------------------------------------------------------
 [ t a b l e . p y ]  |
----------------------+


~$ python3 table.py -i <input_list.txt> -o <output_file.csv> -a <api_key>

--------------------------------------------------------------------------------------------
 [ i n t e r a c t o r s . p y ]  |    ----- table.py output
----------------------------------+    |
                                       |
                                       V
~$ python3 interactors.py -i <input_file.csv> -o <output_file.csv>

--------------------------------------------------------------------------------------------
 [ c o n c a t . p y ]  |          ----- table.py or interactors.py output
------------------------+          |         (two or more allowed as input to concat.py)
                                   |
                                   V
~$ python3 concat.py -i <input_file.csv> <input_file2.csv> -o <output_file.csv>

--------------------------------------------------------------------------------------------
 [ g r a p h . p y ]  |            ----- table.py, interactors.py, or concat.py output
----------------------+            |
                                   |
                                   V
~$ python3 graph.py -i <input_file.csv> -l <labels_option> -m <mode> -o <output_file.png>

graph.py options:

  -h, --help            show this help message and exit

  -m MODE, --mode MODE  argumnents: "gpn" or "protein_interactions"
                        (genotype/phenotype or protein interactions type)

  -i INPUT, --input INPUT
                        <INPUT_FILENAME.csv> (Input table)

  -l LABELS, --labels LABELS
                        arguments: "subtype", "overlapping", "interactors", or "all"

  -o OUTPUT, --output OUTPUT
                        <OUTPUT_FILENAME.png>


GenoPheno (c) 2022-01-27 by Devin Keane
