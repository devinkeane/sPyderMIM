#!/bin/bash

cat << 'EOF'



                           _      _ _ () _ _
                         s||)yder//\/\[]//\/\
                          ||



EOF

sleep 2

python3 split_list.py -i $1 -o $3 > >(tee EXECUTION_OUTPUT.txt) 2>&1

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
  $TABLE_COMMAND > >(tee -a EXECUTION_OUTPUT.txt) 2>&1
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
$CONCAT_COMMAND > >(tee -a EXECUTION_OUTPUT.txt) 2>&1


CONCAT_COMMAND="python3 concat.py -i"
ITERATION=0
while read p; do
  CONCAT_COMMAND+=" ";CONCAT_COMMAND+="$p"; ITERATION=$(($ITERATION + 1))
done<genes_tables_list.txt
CONCAT_COMMAND+=" -o "
CONCAT_COMMAND+=$3_genes_concatenated.csv
$CONCAT_COMMAND > >(tee -a EXECUTION_OUTPUT.txt) 2>&1

python3 interactors.py -i $3_genes_concatenated.csv -m omim -o $3_interactors.csv > >(tee -a EXECUTION_OUTPUT.txt) 2>&1

python3 graph.py -i $3_interactors.csv -m interactors -l all -o $3_interactors > >(tee -a EXECUTION_OUTPUT.txt) 2>&1
python3 graph.py -i $3_interactors.csv -m omim_genes_interactions -g $3_genes_concatenated.csv -l all -o $3_omim_genes_interactions > >(tee -a EXECUTION_OUTPUT.txt) 2>&1
python3 graph.py -i $3_genes_concatenated.csv -m omim_genes -l all -o $3_genes > >(tee -a EXECUTION_OUTPUT.txt) 2>&1
python3 graph.py -i $3_clinical-features_concatenated.csv -m omim_features -l all -o $3_clinical-features > >(tee -a EXECUTION_OUTPUT.txt) 2>&1


python3 enrichment.py -i $3_genes_concatenated.csv -m omim -o $3_genes_enrichment.csv > >(tee -a EXECUTION_OUTPUT.txt) 2>&1
python3 enrichment.py -i $3_interactors.csv  -m interactors -o $3_interactors_enrichment.csv > >(tee -a EXECUTION_OUTPUT.txt) 2>&1

NOW=$( date '+%F_%H%M%S' )
mkdir $3_SPYDERMIM_RESULTS_$NOW

cat EXECUTION_OUTPUT.txt | sed -e 's/.*\r//' >>EXECUTION_OUTPUT_$NOW.txt
rm EXECUTION_OUTPUT.txt

mv $3_clinical-features_concatenated.csv $3_genes_enrichment.csv $3_genes_NETWORK_SUMMARY.txt $3_genes.graphml $3_genes.png $3_genes_concatenated.csv $3_interactors.csv $3_interactors.png $3_omim_genes_interactions.png $3_omim_genes_interactions_NETWORK_SUMMARY.txt $3_omim_genes_interactions.graphml $3_clinical-features.png $3_interactors_enrichment.csv $3_clinical-features_NETWORK_SUMMARY.txt $3_clinical-features.graphml $3_interactors_NETWORK_SUMMARY.txt $3_interactors.graphml EXECUTION_OUTPUT_$NOW.txt $3_SPYDERMIM_RESULTS_$NOW

rm -R ./$3_MIM_directory/
rm -R ./$3_separate_tables/
rm MIM_list.txt
rm genes_tables_list.txt
rm features_tables_list.txt
