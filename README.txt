 ___   ____   __    ___   _      ____
| |_) | |_   / /\  | | \ | |\/| | |_
|_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__
                               ₲Ɇ₦Ø₱ⱧɆ₦Ø v5.3
--------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------
  .-. . . .  . .  . .-. .-. . . |
  `-. | | |\/| |\/| |-| |(   |  |
  `-' `-' '  ` '  ` ` ' ' '  `  |                            
--------------------------------+

GenoPheno is a workflow suite and genotype/phenotype network generator that uses the OMIM
(Online Mendelian Inheritance in Man) database in order to construct a relationship graph that
links diseases by genes and phenotypic outcomes.  The workflow also utilizes the IntAct API,
which allows the user to find the protein products for the genes associated with each MIM number
and the proteins known to interact with each of these.  Additionally, GenoPheno uses the
ToppGene API in order to perform enrichment analysis on all of the genes that encode for
these products.  The user only needs to provide an obtained OMIM API key and a list of OMIM
reference numbers ("phenotypic MIM numbers").

The output features:

   * OMIM genotype/phenotype edge list table (.csv)
   * Intact protein interactors edge list table (.csv)
   * High resolution OMIM genotype/phenotype labeled graph (.png)
   * High resolution Intact protein interactors labeled graph (.png)
   * Genotype/phenotype network statistical analysis (.txt)
   * Protein interactor network statistical analysis (.txt)
   * Genotype/phenotype network graph exchange XML (.gexf)
   * Protein interactor network graph exchange XML (.gexf)
   * Primary gene (OMIM genes) enrichment analysis (.csv)
   * Protein interactor (IntAct genes) enrichment analysis (.csv)

Edge list and .gexf output can be used in other programs.  However, an ultimate goal of our
software is to provide a purely command line based workflow that can allow for graph theory
analysis to be upscaled and applied to larger data sets than could be handled in GUI-based
programs.

_________________________________________________________________________________________________
   __                                           |
  |__)_ _  _  _.|_ _  _    |\/| _ _ |_  _ _ _   |
  | \(-|_)(_)_)||_(_)| \/  |  |(-||||_)(-| _)   |
       |               /                        |
------------------------------------------------+

                                                 ⊂(◉‿◉)つ
 [ s p l i t _ l i s t . p y] will split any MIM list that is greater than 20 into separate
 lists of 20 or less, which can then be used as input for table.py.  This is necessary because
 OMIM API calls are limited to 20 MIM numbers per request.  GenoPheno.sh will perform
 this step automatically.

                                                 (⌐⊙_⊙)
 [ t a b l e . p y ] takes in a .txt file that is a list of "Phenotype MIM numbers," each of
 which correspond to a disease subtype listed in the OMIM Database, collected by the user.  A
 disease subtype in the OMIM database is defined by a gene that is known to give rise to a set
 of phenotypes, which are classified together under the clinical data for that subtype in the
 OMIM database.  The only required input for this workflow is a .txt file that is a simple list
 of phenotype MIM numbers which the user has collected for diseases they wish to investigate
 each separated by a new line (maximum of 20 mims per file, due to API call limits for table.py
 and 5000 total mims when using GenoPheno.sh).

                                               ༼ つ ╹ ╹ ༽つ
 [ i n t e r a c t o r s . p y ]  uses the output from table.py to find the protein products
 of each gene and their protein interactors.

                                                 (✿◠‿◠)
 [ c o n c a t . p y ] will automatically combine multiple .csv outputs from either
 table.py or interactors.py. The output can be used as as single input for graph.py.

                                                ʕっ•ᴥ•ʔっ
 [ g r a p h . p y ] is able to take in the .csv output from either table.py or interactors.py
 and generate a network graph of phenotypes or protein interactors, respectively.  The output is
 a high resolution .png visualization of your network, a .gexf export of your network, and a
 summary .txt file featuring multiple calculations:

        1) top 20 nodes by connectivity degree
        2) top 20 nodes by betweenness centrality
        3) the top 10 nodes in each modularity class by Eigenvector centrality.
        4) network diameter of largest component
        5) triadic closure
        6) number of total nodes and edges        ,,,
                                                \(ʘ‿ʘ)/
 [ e n r i c h m e n t . s h ]  creates a unique list of all interactors found from your
 output table from interactors.py, then performs enrichment analysis on the list using
 the ToppGene enrichment API.  Results have been automatically filtered and sorted in
 ascending order by FDR/B&H Q value < 10e-6.                    ^
                                                                ) )
                                             /\__/\,''`'``'''`'; /
                                            (Ф͡_ᴥ_Ф͡ )_,       ..,;
[ c o n v e r t _ i d s . p y ] will convert any mixed or non-mixed set of genes IDs from
virtually any database identifier type to virtually any other type.  See
"README_conversion_IDs_list.TXT" for an exhaustive list of accepted IDs.


                                             ─=≡Σ((( つ◕ل͜◕)つ

 [ G e n o P h e n o . s h ]  will run the whole workflow automatically on a list of up to
 5,000 MIMs.  Just provide the name of your .txt MIM list file, your OMIM API key, and your
 desired project name.

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
** IMPORTANT NOTES! **   <-- Read first to ensure functionality
----------------------

    ->  GenoPheno.py, the original version of the program, has been replaced by
    GenoPheno.sh.  The whole workflow can now be ran entirely through GenoPheno.sh,
    but you can also use any of the programs in the suite individually or include
    them in your own workflow script!

           ~$ GenoPheno.sh <YOUR_OMIM_API_KEY> <YOUR_PROJECT_NAME>

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

    ->  Do not use phenotypic MIM numbers that correspond to multiple Gene/Locus MIM
    numbers.  Though uncommon, sometimes MIMs are organized this way, and using this
    kind of MIM may cause the program to crash.  You must first access the OMIM entries
    for each of these Gene/Locus MIM numbers individually.  From there, you can select
    and use the phenotype MIM numbers for the subtypes that are listed under these
    entries.

    ->  For MIM numbers where "susceptibility to" is included in the subtype title,
    this program will fetch the ENSEMBL and gene IDs.  However, these MIMs tend to lack
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

(e.g.)
~$ vim input_list.txt

------
615291              [input_list.txt]
614505                      |
120580      <---------------+
and so on ...
------
--------------------------------------------------------------------------------------------

3) Set up an Anaconda environment with the necessary dependencies (Anaconda required):

~$ conda create -n GenoPheno python=3.9 scipy pandas matplotlib curl networkx
~$ conda activate GenoPheno

--------------------------------------------------------------------------------------------

4) Execute with the following syntax/options:

--------------------------------------------------------------------------------------------
 [ G e n o P h e n o . s h ]    |     You can run your entire workflow with one command
--------------------------------+     using GenoPheno.sh with a list of MIMs that is greater
  -. .-.   .-. .-.   .-. .-.   .      than 20.  Alternatively, you can use separate programs
 ||\|||\ /|||\|||\ /|||\|||\ /|       in the suite individually, but table.py will require
 |/ \|||\|||/ \|||\|||/ \|||\||       that your MIM list has 20 or less MIMS, due to the
 ~   `-~ `-`   `-~ `-`   `-~ `-       number of MIMs allowed per API call.  split_lists.py
                                      will create a directory and split your list into multiple
  20 MIM lists if you are using individual programs in the suite.  Read further for usage on
  all programs included in this repository.

                                   +---------- list of OMIM MIM numbers, maximum of 5,000
                                   |
                                   |
Input:                             V

~$ ./GenoPheno.sh <large_input_mim_list.txt> <OMIM_API_key> <project_name>

                                                   A                A
                                                   |                |
        Your API key obtained from OMIM -----------+                |
                                                                    |
        Project name for automatically naming multiple files  ------+

Output:

   * project_name_concatentated.csv  (OMIM genotype/phenotype table)
   * project_name_interactors.csv  (Intact protein interactors table)
   * project_name_concatentated.png  (High resolution OMIM genotype/phenotype labeled graph)
   * project_name_interactors.png  (High resolution Intact protein interactors labeled graph)
   * project_name_phenotypes_NETWORK_SUMMARY.txt  (Genotype/phenotype network statistical analysis)
   * project_name_interactors_NETWORK_SUMMARY.txt  (Protein interactor network statistical analysis)
   * project_name_phenotypes.gexf  (Genotype/phenotype network graph exchange XML)
   * project_name_interactors.gexf  (Protein interactor network graph exchange XML)
   * project_name_primary_genes_enrichment.csv  (OMIM genes enrichment analysis)
   * project_name_interactors_enrichment.csv  (Protein interactors enrichment analysis)

--------------------------------------------------------------------------------------------
 [ s p l i t _ l i s t . p y ]  |   Due to OMIM API call limits, GenoPheno.sh uses
--------------------------------+   split_list.py in order to divide your MIM input list
                                    into separate lists of 20 MIMs or less, which can then
                                    be used as input for table.py.  If you are using the
				  programs in the repository individually for your own
				  workflow or script, you can use split_list.py to divide
				  your MIM list into separate input MIM lists for table.py.


~$ python3 split_list.py -i <big_mim_list.txt> -o <project_name>

Output:
    * ./project_name_MIM_directory/project_name0.txt
    * ./project_name_MIM_directory/project_name1.txt
    * ./project_name_MIM_directory/project_name2.txt
    * etc...
--------------------------------------------------------------------------------------------
 [ t a b l e . p y ]  |
----------------------+


~$ python3 table.py -i <./project_name_MIM_directory/project_name0.txt> -o <output_file.csv> -a <api_key>
~$ python3 table.py -i <./project_name_MIM_directory/project_name1.txt> -o <output_file.csv> -a <api_key>
~$ python3 table.py -i <./project_name_MIM_directoryproject_name2.txt> -o <output_file.csv> -a <api_key>

                                                                                  A
                                                                                  |
-------------------------------------------------------------------------------   |   ------
 [ c o n c a t . p y ]  |          +---- table.py or interactors.py output  ------+
------------------------+          |          (two or more files allowed as input to concat.py)
                                   |
                                   V
~$ python3 concat.py -i <input_file.csv> <input_file2.csv> <and_so_on...> -o <output_file.csv>

--------------------------------------------------------------------------------------------
 [ c o n v e r t _ i d s . p y ]  |      "convert_ids.py" will allows you to convert gene
----------------------------------+      IDs from virtually any type (e.g., ENSG, HGNC, PDB,
                                         etc.) to virtually any type that you desire.  Got a
mixed list of IDs of all different types of IDs?  No problem!  Simply create a .txt file with
no headers and enter each ID, separated by a new line.

You can specify which type of symbol you want to convert your list to by specifying the "-s"
option.  There is no need to specify the type of IDs since this is determined automatically,
allowing for your input list to consist of any mixture of types.  IDs that were not successfully
converted will be reported in an output file.

For a full list of accepted IDs, reference "README_conversion_IDs_list.TXT".


~$ python3 concat.py -i <gene_id_list.txt> -s <SYMBOL_TO_CONVERT_TO> -o <gene_id_list_converted.txt>
                                                       A
Input gene ID list example:                            |
                                                       |
------                                                 +----------  Required!
ENSG00000084674                [gene_id_list.txt]
CLEC2B                                |
CUL2                                  |
ENSG00000185176                       |
ENSG00000117713       <---------------+
1ANI
and so on ...
------

Output:
       *  <OUTPUT_FILENAME.txt>                               <--------  List of converted IDs
       *  <OUTPUT_FILENAME>_<TARGET_SYMBOL>_FAILED_IDs.txt    <--------  IDs that failed to convert

--------------------------------------------------------------------------------------------
 [ i n t e r a c t o r s . p y ]  |    +---- table.py output (.csv) ("omim" mode)
----------------------------------+    |               or
                                       |     any_list.txt (text file of gene names)
                                       V
~$ python3 interactors.py -i <input_file> -m <MODE> -o <output_file.csv>
                                                A
                                                |
                                                |
        "list" or "omim" (required) ____________+

interactors.py takes requires one of two modes, specified using the "-m" option:

    1) "omim"
                  "omim" mode finds protein interactors for all of the genes found in
                  the table that table.py produces for output.

    2) "list"
                  "list" mode takes in a .txt file that you have created.  The list can
                  contain mixed types of IDs from most databases, because interactors.py
                  will use convert_ids.py to detect your ID types and convert them to the
                  necessary format.  Create a list like the one in the example for
                  "convert_ids.py" above.

Output:
       *  For "list" mode:
              *   <OUTPUT_FILENAME>.csv
              *   <OUTPUT_FILENAME>_FAILED_IDS.txt

       *  For "omim" mode:
              *   <OUTPUT_FILENAME>.csv

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


Output:
    * ./project_name_NETWORK_SUMMARY.txt
    * ./project_name.png
    * ./project_name.gexf

--------------------------------------------------------------------------------------------
 [ e n r i c h m e n t . p y ]  |       +---- interactors.py output
--------------------------------+       |
                                        |
                                        V

~$ python3 enrichment.py -i <input_file.csv> -o <output_file.csv>

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

				    (c) 2022-01-27 Devin Keane / Feltus Lab
                                   Department of Genetics & Biochemistry
                                   Clemson University