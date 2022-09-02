     ___   ____   __    ___   _      ____                          __,,,,_
    | |_) | |_   / /\  | | \ | |\/| | |_               _ __..-;''`--/'/ /.',-`-.
    |_| \ |_|__ /_/--\ |_|_/ |_|  | |_|__           (`/' ` |  \ \ \\ / / / / .-'/`,_
                                 ₲Ɇ₦Ø₱ⱧɆ₦Ø v7.0    /'`\ \   |  \ | \| // // / -.,/_,'-,
                                                  /<7' ;  \ \  | ; ||/ /| | \/    |`-/,/-.,_,/')
              C l e m s o n                      /  _.-, `,-\,__|  _-| / \ \/|_/  |    '-/.;.\'
                                                 `-`  f/ ;      / __/ \__ `/ |__/ |
                      U n i v e r s i t y             `-'      |  -| =|\_  \  |-' |
                                                          _ __/   /_..-' `  ),'  //
                                                        fL ((__.-'((___..-'' \__.'
    --------------------------------------------------------------------------------------------

    --------------------------------------------------------------------------------------------
      .-. . . .  . .  . .-. .-. . . |
      `-. | | |\/| |\/| |-| |(   |  |
      `-' `-' '  ` '  ` ` ' ' '  `  |                            
    --------------------------------+

<span style="font-family:Courier">

<font color="lime">GenoPheno is an automated workflow for in silico hypothesis testing.  It also serves as a
"Swiss army knife" of genomics analysis tools.  A user is free to utilize any of the programs
in the repository, but they may also just use GenoPheno.sh to execute most of the suite altogether
as one easily automated workflow and in silico hypothesis testing platform.  The user only needs to
provide API key obtained from OMIM and a list of OMIM reference numbers in a text file ("phenotypic
MIM numbers" beginning with "#" prefix).</font>

<font color="lime">GenoPheno.sh, the full workflow execution program, is a genotype/phenotype network generator that uses
the OMIM (Online Mendelian Inheritance in Man) database in order to construct a relationship graph that
links diseases by genes and phenotypic outcomes.  The workflow also utilizes the IntAct API, which
allows the user to find the protein products for the genes associated with each MIM number and the
proteins \known to interact with each of these.  Additionally, GenoPheno uses the ToppGene API in order
to perform enrichment analysis on all the genes that encode for these products.  By the end of the
workflow, GenoPheno.sh automatically generates two enrichment analyses:</font>

    1)  for the genes obtained from OMIM

    2)  for the genes encoding for all proteins that interact with each of the OMIM Genes


<font size="4"><font color="aqua"> A single workflow run using GenoPheno.sh features the following output, saved to a date/time-stamped folder: </font></font>

     * OMIM gene edge list table (.csv)
     * OMIM clinical features edge list table (.csv)
     * IntAct protein interactors edge list table (.csv)
     * High resolution labeled network image of gene overlap between OMIM disorders (.png)
     * High resolution labeled network image of genes interactions from OMIM disorders (.png)
     * High resolution OMIM clinical data labeled network image (.png)
     * High resolution IntAct protein interactors labeled network image (.png)
     * OMIM gene network statistical analysis (.txt)
     * OMIM clinical features network statistical analysis (.txt)
     * Protein interactor network statistical analysis (.txt)
     * OMIM gene network graph exchange file XML (.gexf)
     * OMIM clinical features network graph exchange file XML (.gexf)
     * Protein interactor network graph exchange file XML (.gexf)
     * Primary gene (OMIM genes) enrichment analysis (.csv)
     * Protein interactor (IntAct genes) enrichment analysis (.csv)


<font color="lime"> Edge list and .gexf output can be used in other programs.  However, an ultimate goal of our
software is to provide a purely command line based workflow that can allow for graph theory
analysis to be upscaled and applied to larger data sets than could be handled in GUI-based
programs. </font>

---

<font size="5"><font color="aqua">The major planned objectives for future versions of GenoPheno include:</font></font>

    *   Adding gene regulatory network discovery capability using GTEX eQTLs

           *   prediction of upstream and downstream chromatin regulatory elements
           *   mapping gene cross-talk directionally
           *   linking DNA sequence with RNA expression across tissues

    *   Adding potential drug discovery functionality




---

    _________________________________________________________________________________________________
       __                                           |
      |__)_ _  _  _.|_ _  _    |\/| _ _ |_  _ _ _   |
      | \(-|_)(_)_)||_(_)| \/  |  |(-||||_)(-| _)   |
           |               /                        |
    ------------------------------------------------+

                                                 ⊂(◉‿◉)つ
--- 
<font size="3"><font color="aqua"><b>[ s p l i t _ l i s t . p y ]</b></span></font></font> <font color="lime">will split any MIM list that is greater than 20 into separate
 lists of 20 or less, which can then be used as input for table.py.  This is necessary because
 OMIM API calls are limited to 20 MIM numbers per request.  GenoPheno.sh will perform
 this step automatically.

---
                                                 (⌐⊙_⊙)
 <font size="3"><font color="aqua"><b>[ t a b l e . p y ]</b></span></font></font> takes in a .txt file that is a list of "Phenotype MIM numbers," each of
 which correspond to a disease subtype listed in the OMIM Database, collected by the user.  A
 disease subtype in the OMIM database is defined by a gene that is known to give rise to a set
 of phenotypes, which are classified together under the clinical data for that subtype in the
 OMIM database.  The only required input for this workflow is a .txt file that is a simple list
 of phenotype MIM numbers which the user has collected for diseases they wish to investigate
 each separated by a new line (maximum of 20 mims per file, due to API call limits for table.py
 and 5000 total mims when using GenoPheno.sh).

table.py uses the [OMIM rest API](https://www.omim.org/help/api) to retreive genetic and clinical data.

---
                                               ༼ つ ╹ ╹ ༽つ
 <font size="3"><font color="aqua"><b>[ i n t e r a c t o r s . p y ]</b></span></font></font>  uses the output from table.py to find the protein products
 of each gene and their protein interactors.


interactors.py uses EBI's [IntAct API](https://www.ebi.ac.uk/intact/documentation/technical_corner)

---
                                                 (✿◠‿◠)
 <font size="3"><font color="aqua"><b>[ c o n c a t . p y ]</b></span></font></font> will automatically combine multiple .csv outputs from either
 table.py or interactors.py. The output can be used as as single input for graph.py.

---
                                                ʕっ•ᴥ•ʔっ
 <font size="3"><font color="aqua"><b>[ g r a p h . p y ]</b></span></font></font> is able to take in the .csv output from either table.py or interactors.py
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

---

graph.py uses the [NetworkX](https://networkx.org/documentation/stable/index.html) Python library to create graph objects.

---

<font size="3"><font color="aqua"><b>[ e n r i c h m e n t . p y ]</b></span></font></font>  creates a unique list of all interactors found from your
 output table from interactors.py, then performs enrichment analysis on the list using
 the ToppGene enrichment API.  Results have been automatically filtered and sorted in
 ascending order by FDR/B&H Q value < 10e-6.

enrichment.py utilizes the Computational Medicine Center's [ToppGene API](https://toppgene.cchmc.org/API/) in order to perform enrichment analysis.
 
                                                                 ^
                                                                ) )
                                             /\__/\,''`'``'''`'; /
                                            (Ф͡_ᴥ_Ф͡ )_,       ..,;

<font size="3"><font color="aqua"><b>[ c o n v e r t _ i d s . p y ]</b></span></font></font> will convert any mixed or non-mixed set of genes IDs from
virtually any database identifier type to virtually any other type.  See
[README_conversion_IDs_list,txt](README_conversion_IDs_list.txt) for an exhaustive list of accepted IDs.

convert_ids.py performs gene name conversions by using g:Profiler's [ g:Convert API](https://biit.cs.ut.ee/gprofiler/page/apis).


                                             ─=≡Σ((( つ◕ل͜◕)つ

<font size="3"><font color="aqua"><b>[ G e n o P h e n o . s h ]</b></span></font></font>  will run the whole workflow automatically on a list of up to
 5,000 MIMs.  Just provide the name of your .txt MIM list file, your OMIM API key, and your
 desired project name. </font>

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
<font size="5"><font color="lime">IMPORTANT NOTES!!!</font></font>

<font size="4"><font color="aqua">** Read first to ensure functionality **</font></font>
___
<font color="lime">* GenoPheno.sh is the main program that will execute the entire workflow automatically.
Additionally, you can also use any of the programs in the suite individually or include
them in your own workflow script!</font>


       ~$ GenoPheno.sh <MIM_LIST.TXT> <YOUR_OMIM_API_KEY> <YOUR_PROJECT_NAME>  (do not use file
                                                                              extension for
                                                                              project name)


<font color="lime">* GenoPheno is designed to take in lists of "Phenotype MIM numbers" as the starting input of the workflow.
DO NOT use "Gene/Locus MIM numbers" as these will cause the program to crash.  This program does not process
MIM numbers with the "+", "%", "*", or "^" prefixes as these are irrelevent to the objectives of this software
and the program will crash if these types of MIM numbers are used.</font>




<font color="lime">Following the instructions above will ensure optimal program functionality and output.</font>

--------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
      |   |
      |   |   __|   _` |   _` |   _ \
      |   | \__ \  (   |  (   |   __/
      \___/ ____/ \__,_| \__, | \___|
                         |___/
--------------------------------------------------------------------------------------------

1) <font size="3"><font color="lime">Obtain an API key through OMIM:  https://www.omim.org/api </font></font>

--------------------------------------------------------------------------------------------

2) <font size="3"><font color="lime">Create a list of OMIM reference ids ("MIM" numbers)

(20 MAXIMUM if input is for table.py, 5000 MAXIMUM if using GenoPheno.sh): </font></font>


~$  ` vim input_list.txt `
    
    ------
    615291              [input_list.txt]
    614505                      |
    120580      <---------------+
    and so on ...
    ------
--------------------------------------------------------------------------------------------

3) <font size="3"><font color="lime">Set up an Anaconda environment with the necessary dependencies (Anaconda required): </font></font>


~$ ` conda create -n GenoPheno python=3.9 scipy=1.7 pandas matplotlib curl networkx requests`

~$ ` conda activate GenoPheno` 

--------------------------------------------------------------------------------------------

4) <font size="3"><font color="lime">Execute with the following syntax/options: </font></font>
    
--------------------------------------------------------------------------------------------

    
     --------------------------------------------------------------------------------------------
     [ G e n o P h e n o . s h ]    |     You can run your entire workflow with one command
     -------------------------------+     using GenoPheno.sh with a list of MIMs that is greater
     -. .-.   .-. .-.   .-. .-.   .       than 20.  Alternatively, you can use separate programs
     ||\|||\ /|||\|||\ /|||\|||\ /|       in the suite individually, but table.py will require
     |/ \|||\|||/ \|||\|||/ \|||\||       that your MIM list has 20 or less MIMS, due to the
     ~   '-~ '-'   '-~ '-'   '-~ '-       number of MIMs allowed per API call.  split_lists.py
                                          will create a directory and split your list into multiple
      20 MIM lists if you are using individual programs in the suite.  Read further for usage on
      all programs included in this repository.

                                   +---------- list of OMIM MIM numbers, maximum of 5,000
                                   |
                                   |
    Input:                         V

    ~$  ./GenoPheno.sh <large_input_mim_list.txt> <OMIM_API_key> <project_name>

                                                       A                A
                                                       |                |
            Your API key obtained from OMIM -----------+                |
                                                                        |
            Project name for automatically naming multiple files  ------+


<font color="lime">Output:  (all results are moved to a new folder
          named after the project name, followed the a date/timestamp)

   * project_name_genes.csv  (OMIM geno table)
   * project_name_clinical-features.csv  (OMIM clinical features table)
   * project_name_interactors.csv  (IntAct protein interactors table)
   * project_name_genes.png  (High resolution OMIM gene labeled graph)
   * project_name_omim_genes_interactions.png  (High resolution labeled graph of OMIM disease gene interactions)
   * project_name_clinical-features.png  (High resolution OMIM clinical data labeled graph)
   * project_name_interactors.png  (High resolution IntAct protein interactors labeled graph)
   * project_name_clinical-features_NETWORK_SUMMARY.txt  (OMIM clinical features network statistical analysis)
   * project_name_genes_NETWORK_SUMMARY.txt  (OMIM genes network statistical analysis)
   * project_name_omim_genes_interactions_NETWORK_SUMMARY.txt  (OMIM gene interactions network statistical analysis)
   * project_name_interactors_NETWORK_SUMMARY.txt  (IntAct protein interactor network statistical analysis)
   * project_name_clinical-features.gexf  (OMIM clinical features network graph exchange XML)
   * project_name_genes.gexf  (OMIM disease gene overlap network graph exchange XML)
   * project_name_omim_genes_interactions.gexf  (OMIM disease gene interactions network graph exchange XML)
   * project_name_interactors.gexf  (Protein interactor network graph exchange XML)
   * project_name_genes_enrichment.csv  (OMIM genes enrichment analysis)
   * project_name_interactors_enrichment.csv  (IntAct protein interactors enrichment analysis) </font>

    --------------------------------------------------------------------------------------------
     [ s p l i t _ l i s t . p y ]  |   Due to OMIM API call limits, GenoPheno.sh uses
    --------------------------------+   split_list.py in order to divide your MIM input list
                                        into separate lists of 20 MIMs or less, which can then
                                        be used as input for table.py.  If you are using the
			    	  programs in the repository individually for your own
			    	  workflow or script, you can use split_list.py to divide
			    	  your MIM list into separate input MIM lists for table.py.


    ~$ python3 split_list.py -i <big_mim_list.txt> -o <project_name>   <-- (no file extension needed)

    Output:
    * ./project_name_MIM_directory/project_name0.txt
    * ./project_name_MIM_directory/project_name1.txt ---+
    * ./project_name_MIM_directory/project_name2.txt    |
    * etc...                                            |
                                                        |
                                                        |
    -------------------------------------------------   |   -----------------------------------
     [ t a b l e . p y ]  |                             |
    ----------------------+                             |
                                                        V

    python3 table.py -i <./project_name_MIM_directory/project_name0.txt> -o <output_name_0> -a <api_key>
    python3 table.py -i <./project_name_MIM_directory/project_name1.txt> -o <output_name_1> -a <api_key>
    python3 table.py -i <./project_name_MIM_directory/project_name2.txt> -o <output_name_2> -a <api_key>
    
        Output: |     * output_name_genes_0.csv         * output_name_clinical-features_0.csv                                         
        --------+     * output_name_genes_1.csv         * output_name_clinical-features_1.csv
                      * output_name_genes_2.csv         * output_name_clinical-features_2.csv
                                                                                      A
                                                                                      |
    -------------------------------------------------------------------------------   |   ------
     [ c o n c a t . p y ]  |          +---- table.py or interactors.py output  ------+
    ------------------------+          |          (two or more files allowed as input to concat.py)
                                       |
                                       V
    ~$ python3 concat.py -i <input_file_0.csv> <input_file_1.csv> <and_so_on...> -o <output_file.csv>`


    --------------------------------------------------------------------------------------------
     [ i n t e r a c t o r s . p y ]  |    +---- table.py output (.csv) ("omim" mode)
    ----------------------------------+    |               or
                                           |     any_list.txt (text file of gene names)
                                           V
    ~$ python3 interactors.py -i <input_file> -m <MODE> -o <output_name>    <-- (no file extension)
                                                    A
                                                    |
                                                    |
            "list" or "omim" (required) ____________+
    
<font size="3"><font color="lime">interactors.py requires one of two modes, specified using the "-m" option:

### 1) "omim"
* "omim" mode finds protein interactors for all of the genes found in the table that table.py produces for output.

### 2) "list"
* "list" mode takes in a .txt file that you have created.  The list can
    contain mixed types of IDs from most databases, because interactors.py
    will use convert_ids.py to detect your ID types and convert them to the
    necessary format.  Create a list like the one in the example for
     "convert_ids.py" above. </font></font>

<!-- end the list -->

         (e.g.) vim list.txt:
         
         ------
         ENSG00000185176
         CUL2
         1ANI
         and so on...
         ------
    

<font size="3"><font color="lime">Output:
* For "list" mode:
  * <font color="aqua"><OUTPUT_FILENAME>.csv
  * <OUTPUT_FILENAME>_FAILED_IDS.txt</font>


* For "omim" mode: </font>
  * <font color="aqua"><OUTPUT_FILENAME>.csv </font> </font>

<!-- end the list -->

    --------------------------------------------------------------------------------------------
     [ g r a p h . p y ]  |            +---- table.py, interactors.py, or concat.py output
    ----------------------+            |
                                       |
                                       V
    ~$ python3 graph.py -i <input_file.csv> -m <mode> -l <labels_option> -o <output_file.png>
    
    graph.py options:
    
      -h, --help            show this help message and exit
    
      -m MODE, --mode MODE  "list" (find protein interactors for a .txt file of most gene ID types)
                                                                                                 A
	    				                    or         (see README_conversion_IDs_list.txt)   ___|           
    
	    	                "omim" (find protein interactors for your table.py OMIM genes output edge list)
    
      -i INPUT, --input INPUT

                        <INPUT_FILENAME.csv> (Input table)
    
      -l LABELS, --labels LABELS
    
                            arguments: arguments: "all" or "none"
    
      -o OUTPUT, --output OUTPUT
    
                            <OUTPUT_FILENAME.png>
    
    
<font size="3"><font color="lime">Output:
* ./project_name_NETWORK_SUMMARY.txt
* ./project_name.png
* ./project_name.gexf </font></font>

<!-- end the list -->

    --------------------------------------------------------------------------------------------
     [ e n r i c h m e n t . p y ]  |       +---- interactors.py output
    --------------------------------+       |
                                            |
                                            V
    
    ~$ python3 enrichment.py -i <input_file.csv> -o <output_file.csv>
    
    --------------------------------------------------------------
<font size="3"><font color="lime">Output:
* Enrichment analysis (.csv) </font>
  * <font color="aqua">filtered by FDR/B&H Q value < 10e-6 in ascending order. </font></font>

---

    --------------------------------------------------------------------------------------------
     [ c o n v e r t _ i d s . p y ]  |      "convert_ids.py" will allow you to convert gene
    ----------------------------------+      IDs from virtually any type (e.g., ENSG, HGNC, PDB,
                                         etc.) to virtually any type that you desire.  Got a
    mixed list of IDs of all different types?  No problem!  Simply create a .txt file with one
    column of gene IDs and no header.
    
    You can specify which type of symbol you want to convert your list to by specifying the "-s"
    option.  There is no need to specify the type of IDs since this is determined automatically,
    allowing for your input list to consist of any mixture of types.  IDs that were not successfully
    converted will be reported in an output file.
    
<font size="3"><font color="lime">For a full list of accepted IDs, reference [README_conversion_IDs_list,txt](README_conversion_IDs_list.txt). </font></font>

<font size="3"><font color="aqua"><b>Conversions were made using API calls to the g:Profiler web service: https://biit.cs.ut.ee/gprofiler/gost </font></font></b>

    
    ~$ python3 convert_ids.py -i <gene_id_list.txt> -s <SYMBOL_TO_CONVERT_TO> -o <gene_id_list_converted.txt>
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
<font size="3"><font color="lime">
Output:
* <OUTPUT_FILENAME.txt>
  * <font color="aqua">(List of converted IDs) </font>


* <OUTPUT_FILENAME>_<TARGET_SYMBOL>_FAILED_IDs.txt  
    * <font color="aqua"> (IDs that failed to convert) </font>

</font></font>
<!-- end the list -->

---

    --------------------------------------------------------------------------------------------
     [ c o n t r o l _ t r i a l s _ g e n e r a t o r . p y ]  |    * BETA *
    ------------------------------------------------------------+ 
                                         
    The random controls trials generator was developed to test the ability of GenoPheno
    to find overlap and clustering between genes and clinical data for diseases
    that are thought to be related.
    
    It allows a user working on a supercomputing cluster to specify a set of MIMs of
    interest which they can run through GenoPheno against a specified number of MIMs
    randomly selected from the OMIM database.

    For example, suppose a user runs an experiment with 3 MIMs corresponding to varieties
    of hypophosphatasia against 45 MIMs corresponding to Ehlers-Danlos Syndrome and
    Osteogenesis Imperfectica.  Using the random controls trials generator, the user
    is able to run dozens of experiments testing the three hypophosphatasia MIMs against
    45 randomly selected MIMs.

    Alternatively, the user could specify that they do not want to include any experimental
    MIMs in each batch by pressing "RETURN" when prompted to enter experimental MIMs
    (hypophosphatasia, in the above example).  In this case, the program will generate
    experiments where each batch of MIMs simply contains 48 randomly selected MIMs, rather
    containing set of specified experimental MIMs together with the 45 random MIMs.

    ---------------------------------------------------------------------------------
    (e.g.)
                                  ____
                                 /\' .\    _____
         C o n t r o l          /: \___\  / .  /\
          T r i a l s           \' / . / /____/..\
       G e n e r a t o r         \/___/  \'  '\  /
               ( beta )                   \'__'\/


     (✿◠‿◠)  Please enter a name for your experiment:  TEST_1
     (✿◠‿◠)  Please tell me how many MIM numbers would you like to generate per batch:  48
     (✿◠‿◠)  Please tell me how many batches you would like to generate:  10 
     (✿◠‿◠)  Please list the experimental MIMs to run against controls, separated by spaces
     (or press RETURN for none):  146300 241510 241500
     (✿◠‿◠)  Please enter your API key so I can run the experiment:

   ---------------------------------------------------------------------------------

    In the given example above, 10 .pbs scripts will be generated, each running a GenoPheno
    experiment containing 3 specified MIMs (Hypophosphatasia) and 45 randomly selected MIMs.
    For now, computer cluster parameters cannot be specified.  By default, the .pbs
    scripts will request "fdr" interconnect, 100gb RAM, 20 CPUs, and 60 hours of walltime.

---

<font color="gold">(c) 2022-01-27 Devin Keane / Feltus Lab</font>

<font color="gold">Department of Genetics & Biochemistry</font>

<font color="gold">Clemson University</font>

</span>