 ___   ____   __    ___   _      ____
| |_) | |_   / /\  | | \ | |\/| | |_
|_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__
                               ₲Ɇ₦Ø₱ⱧɆ₦Ø

GenoPheno is a workflow suite and genotype/phenotype network generator that uses the OMIM
(Online Mendelian Inheritance in Man) database in order to construct a relationship graph that
links diseases by genes and phenotypes.  The workflow also utilizes the IntAct API, which
allows the user to find the protein products for the genes associated with each MIM number
and the proteins known to interact with each of these.  Large tables and network graphs can
be constructed in order to investigate large batches of disease subtypes.

table.py takes in a .txt file that is a list of "Phenotype MIM numbers," each of which correspond
to a disease subtype listed in the OMIM Database, collected by the user.  A disease subtype in
the OMIM database is defined by a gene that is known to give rise to a set of phenotypes, which
are classified together under the clinical data for that subtype in the OMIM database.  The only
required input for this workflow is a .txt file(s) that is a simple list of phenotype MIM numbers
which the user has collected for diseases they wish to investigate, each separated by a new line
(maximum of 20 mims per file, due to API call limits, and 5000 total mims per day).

interactors.py uses the output from table.py to find the protein products of each gene and
their protein interactors.

graph.py is able to take in the .csv output from either table.py or interactors.py and
generate a network graph of phenototypes or protein interactors, respectively.

Output from these tables can be used as edge list input for other programs, such as Cytoscape!
However, the goal of this program is to provide a purely command line based interface that can
allow for graph theory analysis to be upscaled and applied to larger data sets than could be
handled in GUI-based programs.

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
** IMPORTANT NOTES FOR VERSION 3.1 ! **   <-- Read first to ensure functionality
---------------------------------------
    ->  GenoPheno.py, the original version of the program, has been replaced by
    table.py, but it will still remain on GitHub for the time being.  table.py
    may later be renamed to GenoPheno.py after the original version is removed.

    ->  GenoPheno only takes in lists of "Phenotype MIM numbers" as the
    starting input of the workflow.  DO NOT use "Gene/Locus MIM numbers" as
    these will cause the program to crash.  Future versions of the program will
    aim to ignore these MIM numbers.

    ->  When creating a .txt list of mim numbers, ONLY USE MIM numbers with the "#"
    prefix.  This program does not process MIM numbers with the "+", "%", or "^"
    prefixes as these are irrelevent to the objectives of this software and the
    program will crash if these types of MIM numbers are used.

    ->  Do not use mim numbers for listings where "susceptibility to" is
    is included in the subtype title.  This program does not support MIM numbers
    for this nor MIM numbers that correspond to individual genes rather than
    subtypes.

    ->  Do not use mim numbers that correspond to multiple genes.  These MIM numbers
    represent families of sub-phenotypes that are each associated with a single gene.
    You must first find the Gene/Locus MIM number associated with this MIM number, then
    select and use the phenotypic MIM numbers for the subtypes that are listed under
    it.

** table.py will likely continue to process your MIM list and problematic MIMs will
be listed in your output.  They should still be avoided if possible.  **
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

3) Set up an Anaconda environment with the necessary dependencies (Anaconda required):

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib networkx
~$ conda activate GenoPheno

** Note: these dependencies reflect an earlier version of the program, so you may need to
install others in your Anaconda environment, depending on your error output.  I am currently
working on updating this for the README file. -- DK 2022-03-03 **

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
~$ python3 graph.py -i <input_file.csv> -m <mode> -l <labels_option> -o <output_file.png>

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


GenoPheno (c) 2022-01-27 by Devin Keane / Clemson University
