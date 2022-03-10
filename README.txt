 ___   ____   __    ___   _      ____
| |_) | |_   / /\  | | \ | |\/| | |_
|_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__
                               ₲Ɇ₦Ø₱ⱧɆ₦Ø 4.0

GenoPheno is a workflow suite and genotype/phenotype network generator that uses the OMIM
(Online Mendelian Inheritance in Man) database in order to construct a relationship graph that
links diseases by genes and phenotypic outcomes.  The workflow also utilizes the IntAct API,
which allows the user to find the protein products for the genes associated with each MIM number
and the proteins known to interact with each of these.  Large tables and network graphs can
be constructed in order to investigate large batches of disease subtypes.  The user only needs
to provide an obtained OMIM API key and a list of OMIM reference numbers ("phenotypic MIM numbers").

These output tables can be used as edge list input for other programs, such as Cytoscape.
However, the goal of our program is to provide a purely command line based interface that can
allow for graph theory analysis to be upscaled and applied to larger data sets than could be
handled in GUI-based programs.

_________________________________________________________________________________________________
                           
_|_  _. |_  |  _    ._     
 |_ (_| |_) | (/_ o |_) \/ takes in a .txt file that is a list of "Phenotype MIM numbers," each of 
                    |   /  which correspond to a disease subtype listed in the OMIM Database,
collected by the user.  A disease subtype in the OMIM database is defined by a gene that is known to give rise to a set of phenotypes, which are classified together under the clinical data for that subtype in the OMIM database.  The only required input for this workflow is a .txt file that is
a simple list of phenotype MIM numbers which the user has collected for diseases they wish to
investigate, each separated by a new line (maximum of 20 mims per file, due to API call limits
for table.py and 5000 total mims when using GenoPheno.sh).

 _ _ |._|_  |. __|_  _   
_\|_)|| | __||_\ | .|_)\/  will split any MIM list that is greater than 20 into separate lists 
  |                 |  /   of 20 or less, which can then be used as input for table.py.
			  GenoPheno.sh will perform this step automatically.       
                      
 o ._ _|_  _  ._ _.  _ _|_  _  ._  _   ._     
 | | | |_ (/_ | (_| (_  |_ (_) |  _> o |_) \/ 
                                       |   /  uses the output from table.py to find the protein
					    products of each gene and their protein interactors.                   
  _  _  ._   _  _. _|_   ._     
 (_ (_) | | (_ (_|  |_ o |_) \/  will automatically combine multiple .csv outputs from either
                         |   /   table.py or interactors.py. The output can be used as as single
			        input for graph.py.
  _  ._ _. ._  |_    ._     
 (_| | (_| |_) | | o |_) \/  is able to take in the .csv output from either table.py or interactors.py 
  _|       |         |   /   and generate a network graph of phenotypes or protein interactors,
			    respectively.
  __              _                           
 /__  _  ._   _  |_) |_   _  ._   _     _ |_  
 \_| (/_ | | (_) |   | | (/_ | | (_) o _> | |  will run the whole workflow automatically on a list of
					     up to 5,000 MIMs, saving one large

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
** IMPORTANT NOTES FOR VERSION 4.0 ! **   <-- Read first to ensure functionality
---------------------------------------

    ->  GenoPheno.py, the original version of the program, has been replaced by
    GenoPheno.sh.  The whole workflow can now be ran entirely through GenoPheno.sh,
    but you can also use any of the programs in the suite individually or to include
    in your own workflow script.

    ->  GenoPheno only takes in lists of "Phenotype MIM numbers" as the
    starting input of the workflow.  DO NOT use "Gene/Locus MIM numbers" as
    these will cause the program to crash.  Even though these MIMs cannot be used, you
    can usually locate phenotype MIMs from the phenotypic series that these MIMs are
    associated with.  See "Usage" below for instructions on how to use GenoPheno.sh as
    well as how to use each individual program in the suite.

    ->  When creating a .txt list of MIM numbers, ONLY USE MIM numbers with the "#"
    prefix.  This program does not process MIM numbers with the "+", "%", "*", or "^"
    prefixes as these are irrelevent to the objectives of this software and the
    program will crash if these types of MIM numbers are used.

    ->  Do not use MIM numbers that correspond to multiple genes.  These MIM numbers
    represent families of sub-phenotypes that are each associated with a single gene.
    You must first find the Gene/Locus MIM number associated with this kind of MIM
    number, then select and use the phenotype MIM numbers for the subtypes that are
    listed under it.

    ->  For MIM numbers where "susceptibility to" is included in the subtype title,
    this program will fetch the ENSML and gene IDs.  However, these MIMs tend to lack
    phenotypic data since they describe genes that indirectly influence individuals'
    susceptibility to a disease, not the genes that give rise to the phenotypes
    themselves.



GenoPheno.sh and table.py should continue to process your MIM list and problematic
MIMs will be listed in your output and ignored if unusable.  Following the instructions
above to avoid problematic MIMs will ensure optimal program functionality and output.


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

2) Create list of OMIM reference ids ("MIM" numbers) (20 MAXIMUM if input is for table.py,
5000 MAXIMUM if using GenoPheno.sh):

(i.e.)
~$ vim input_list.txt
------
615291              [input_list.txt]
614505                      |
120580      <---------------+
and so on ...
------
--------------------------------------------------------------------------------------------

3) Set up an Anaconda environment with the necessary dependencies (Anaconda required):

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib networkx
~$ conda activate GenoPheno

--------------------------------------------------------------------------------------------

4) Execute with the following syntax/options:

--------------------------------------------------------------------------------------------
 [ G e n o P h e n o . s h ]    |     You can run your entire workflow with one command
--------------------------------+     using GenoPheno.sh with a list of MIMs that is greater
                                      than 20.  Alternatively, you can use separate programs
                                      in the suite individually, but table.py will require
                                      that your MIM list has 20 or less MIMS, due to the
				    number of MIMs allowed per API call.  split_lists.py
				    will create a directory and split your list into
				    multiple 20 MIM lists if you are using individual
				    programs in the suite.  Read further for usage on all
				    programs included in this repository.


                                   +---------- list of OMIM MIM numbers, maximum of 5,000
                                   |
                                   V

~$ ./GenoPheno.sh <large_input_mim_list.txt> <OMIM_API_key> <project_name>

                                                   A                A
                                                   |                |
        Your API key obtained from OMIM -----------+                |
                                                                    |
        Project name for automatically naming multiple files  ------+

--------------------------------------------------------------------------------------------
 [ s p l i t _ l i s t . p y ]  |
--------------------------------+

~$ python3 split_list.py -i <big_mim_list.txt> -o <project_name>

Output:
    * ./project_name_separate_lists/project_name0.txt
    * ./project_name_separate_lists/project_name1.txt
    * ./project_name_separate_lists/project_name2.txt
    * etc...
--------------------------------------------------------------------------------------------
 [ t a b l e . p y ]  |
----------------------+


~$ python3 table.py -i <./project_name_separate_lists/project_name0.txt> -o <output_file.csv> -a <api_key>
~$ python3 table.py -i <./project_name_separate_lists/project_name1.txt> -o <output_file.csv> -a <api_key>
~$ python3 table.py -i <./project_name_separate_lists/project_name2.txt> -o <output_file.csv> -a <api_key>

                                                                                  A
                                                                                  |
-------------------------------------------------------------------------------   |   ------
 [ c o n c a t . p y ]  |          +---- table.py or interactors.py output  ------+
------------------------+          |          (two or more files allowed as input to concat.py)
                                   |
                                   V
~$ python3 concat.py -i <input_file.csv> <input_file2.csv> <and_so_on...> -o <output_file.csv>

--------------------------------------------------------------------------------------------
 [ i n t e r a c t o r s . p y ]  |    +---- table.py output
----------------------------------+    |
                                       |
                                       V
~$ python3 interactors.py -i <input_file.csv> -o <output_file.csv>


--------------------------------------------------------------------------------------------
 [ g r a p h . p y ]  |            +---- table.py, interactors.py, or concat.py output
----------------------+            |
                                   |
                                   V
~$ python3 graph.py -i <input_file.csv> -m <mode> -l <labels_option> -o <output_file.png>

graph.py options:

  -h, --help            show this help message and exit

  -m MODE, --mode MODE  "geno" (find protein interactor overlap)

					or 

		       "pheno" (find OMIM phenotypic overlap)

  -i INPUT, --input INPUT

                        <INPUT_FILENAME.csv> (Input table)

  -l LABELS, --labels LABELS

                        arguments: "subtype", "overlapping", "interactors", or "all"

  -o OUTPUT, --output OUTPUT

                        <OUTPUT_FILENAME.png>

--------------------------------------------------------------

                         __,,,,_
          _ __..-;''`--/'/ /.',-`-.
      (`/' ` |  \ \ \\ / / / / .-'/`,_
     /'`\ \   |  \ | \| // // / -.,/_,'-,
    /<7' ;  \ \  | ; ||/ /| | \/    |`-/,/-.,_,/')
   /  _.-, `,-\,__|  _-| / \ \/|_/  |    '-/.;.\'
   `-`  f/ ;      / __/ \__ `/ |__/ |
        `-'      |  -| =|\_  \  |-' |
            _ __/   /_..-' `  ),'  //
          fL ((__.-'((___..-'' \__.'

				    (c) 2022-01-27 Devin Keane / Clemson University