#!/bin/bash

python3 split_list.py -i $1 -o $3 2>&1 | tee >(col -pbf >EXECUTION_OUTPUT_temp.txt)

ls ./$3_MIM_directory > file_list.txt
TABLE_COMMAND="python3 table.py -i "

ITERATION=0

mkdir $3_separate_tables
while read p; do
  TABLE_COMMAND="python3 table.py -i "
  TABLE_COMMAND+="./$3_MIM_directory/$p -a $2 -o ./$3_separate_tables/$3$ITERATION.csv"; ITERATION=$(($ITERATION + 1))
  $TABLE_COMMAND 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
done<file_list.txt

ls ./$3_separate_tables > file_list2.txt

CONCAT_COMMAND="python3 concat.py -i"
ITERATION=0
while read p; do
  CONCAT_COMMAND+=" ";CONCAT_COMMAND+="./$3_separate_tables/$p"; ITERATION=$(($ITERATION + 1))
done<file_list2.txt
CONCAT_COMMAND+=" -o "
CONCAT_COMMAND+=$3_concatenated.csv
$CONCAT_COMMAND 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

python3 interactors.py -i $3_concatenated.csv -m omim -o $3_interactors.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

python3 graph.py -i $3_interactors.csv -m geno -l all -o $3_interactors.png 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 graph.py -i $3_concatenated.csv -m pheno -l all -o $3_phenotypes.png 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

python3 enrichment.py -i $3_concatenated.csv -m primary -o $3_primary_genes_enrichment.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 enrichment.py -i $3_interactors.csv  -m interactors -o $3_interactors_enrichment.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

NOW=$( date '+%F_%H%M%S' )
mkdir $3_GENOPHENO_RESULTS_$NOW

sed -e 's/⌐|■ENSEMBLoID:n/ | ENSEMBL ID: /g' -e 's/( °_ʖ °)_\/¯lCalculation complete!/( °_ʖ °)_\/¯ Calculation complete!/g' -e 's/.*extracted.*/ENSEMBL IDs extracted	  ✔/g' EXECUTION_OUTPUT_temp.txt > EXECUTION_OUTPUT_$NOW.txt
rm EXECUTION_OUTPUT_temp.txt

mv $3_concatenated.csv $3_interactors.csv $3_interactors.png $3_phenotypes.png $3_primary_genes_enrichment.csv $3_interactors_enrichment.csv $3_phenotypes_NETWORK_SUMMARY.txt $3_phenotypes.gexf $3_interactors_NETWORK_SUMMARY.txt $3_interactors.gexf EXECUTION_OUTPUT_$NOW.txt $3_GENOPHENO_RESULTS_$NOW

rm -R ./$3_MIM_directory/
rm -R ./$3_separate_tables/
rm file_list.txt
rm file_list2.txt
