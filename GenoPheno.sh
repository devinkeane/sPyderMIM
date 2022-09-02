#!/bin/bash

python3 split_list.py -i $1 -o $3 2>&1 | tee >(col -pbf >EXECUTION_OUTPUT_temp.txt)

ls ./$3_MIM_directory > MIM_list.txt
TABLE_COMMAND="python3 table.py -i "

ITERATION=0
if [[ -d $3_separate_tables ]]
then
    rm -R $3_separate_tables
fi
mkdir $3_separate_tables

while read p; do
  TABLE_COMMAND="python3 table.py -i "
  TABLE_COMMAND+="./$3_MIM_directory/$p -a $2 -o ./$3_separate_tables/$3_$ITERATION"; ITERATION=$(($ITERATION + 1))
  $TABLE_COMMAND 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
done<MIM_list.txt

ls ./$3_separate_tables/*_clinical-features.csv > features_tables_list.txt
ls ./$3_separate_tables/*_genes.csv > genes_tables_list.txt

CONCAT_COMMAND="python3 concat.py -i"
ITERATION=0
while read p; do
  CONCAT_COMMAND+=" ";CONCAT_COMMAND+="$p"; ITERATION=$(($ITERATION + 1))
done<features_tables_list.txt
CONCAT_COMMAND+=" -o "
CONCAT_COMMAND+=$3_clinical-features_concatenated.csv
$CONCAT_COMMAND 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)


CONCAT_COMMAND="python3 concat.py -i"
ITERATION=0
while read p; do
  CONCAT_COMMAND+=" ";CONCAT_COMMAND+="$p"; ITERATION=$(($ITERATION + 1))
done<genes_tables_list.txt
CONCAT_COMMAND+=" -o "
CONCAT_COMMAND+=$3_genes_concatenated.csv
$CONCAT_COMMAND 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

python3 interactors.py -i $3_genes_concatenated.csv -m omim -o $3_interactors.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

python3 graph.py -i $3_interactors.csv -m interactors -l all -o $3_interactors 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 graph.py -i $3_interactors.csv -m omim_genes_interactions -g $3_genes_concatenated.csv -l all -o $3_omim_genes_interactions 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 graph.py -i $3_genes_concatenated.csv -m omim_genes -l all -o $3_genes 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 graph.py -i $3_clinical-features_concatenated.csv -m omim_features -l all -o $3_clinical-features 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)


python3 enrichment.py -i $3_genes_concatenated.csv -m omim -o $3_genes_enrichment.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)
python3 enrichment.py -i $3_interactors.csv  -m interactors -o $3_interactors_enrichment.csv 2>&1 | tee -a >(col -pbf >>EXECUTION_OUTPUT_temp.txt)

NOW=$( date '+%F_%H%M%S' )
mkdir $3_GENOPHENO_RESULTS_$NOW

sed -e 's/⌐|■ENSEMBLoID:n/ | ENSEMBL ID: /g' -e 's/( °_ʖ °)_\/¯lCalculation complete!/( °_ʖ °)_\/¯ Calculation complete!/g' -e 's/.*extracted.*/ENSEMBL IDs extracted	  ✔/g' EXECUTION_OUTPUT_temp.txt > EXECUTION_OUTPUT_$NOW.txt
rm EXECUTION_OUTPUT_temp.txt

mv $3_clinical-features_concatenated.csv $3_genes_enrichment.csv $3_genes_NETWORK_SUMMARY.txt $3_genes.gexf $3_genes.png $3_genes_concatenated.csv $3_interactors.csv $3_interactors.png $3_omim_genes_interactions.png $3_omim_genes_interactions_NETWORK_SUMMARY.txt $3_omim_genes_interactions.gexf $3_clinical-features.png $3_interactors_enrichment.csv $3_clinical-features_NETWORK_SUMMARY.txt $3_clinical-features.gexf $3_interactors_NETWORK_SUMMARY.txt $3_interactors.gexf EXECUTION_OUTPUT_$NOW.txt $3_GENOPHENO_RESULTS_$NOW

rm -R ./$3_MIM_directory/
rm -R ./$3_separate_tables/
rm MIM_list.txt
rm genes_tables_list.txt
rm features_tables_list.txt
