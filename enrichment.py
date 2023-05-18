#                              d8,        d8b
#                             `8P         ?88                                        d8P
#                                          88b                                    d888888P |
#   d8888b  88bd88b   88bd88b  88b d8888b  888888b   88bd8b,d88b  d8888b  88bd88b   ?88'   |
#  d8b_,dP  88P' ?8b  88P'  `  88Pd8P' `P  88P `?8b  88P'`?8P'?8bd8b_,dP  88P' ?8b  88P    |
#  88b     d88   88P d88      d88 88b     d88   88P d88  d88  88P88b     d88   88P  88b    |
#  `?888P'd88'   88bd88'     d88' `?888P'd88'   88bd88' d88'  88b`?888P'd88'   88b  `?8b   |
#                                                                                          |
#  ----------------------------------------------------------------------------------------+
#                                                        丂尸丫ᗪ🝗尺爪工爪  7.3.2
#
#                                                        Created: 2022-03-15
#
#  -----------------------------------------------------------------------------------------

import os
import re
import json
import requests
import csv
import time
import shutil



logo = """
                              d8,        d8b
                             `8P         ?88                                        d8P
                                          88b                                    d888888P |
   d8888b  88bd88b   88bd88b  88b d8888b  888888b   88bd8b,d88b  d8888b  88bd88b   ?88'   |
  d8b_,dP  88P' ?8b  88P'  `  88Pd8P' `P  88P `?8b  88P'`?8P'?8bd8b_,dP  88P' ?8b  88P    |
  88b     d88   88P d88      d88 88b     d88   88P d88  d88  88P88b     d88   88P  88b    |
  `?888P'd88'   88bd88'     d88' `?888P'd88'   88bd88' d88'  88b`?888P'd88'   88b  `?8b   |
                                                                                          |
  ----------------------------------------------------------------------------------------+
                                                        丂尸丫ᗪ🝗尺爪工爪  v7.4.0
"""

print(logo)

OUTPUT_DIR = 'enrichment_results'

def process_all_files(gene_set_library):
    filenames = [filename for filename in os.listdir('.') if re.match(r'.+(genes|interactors|omim_genes_interactions)_NETWORK_SUMMARY\.txt$', filename)]
    for filename in filenames:
        process_file(gene_set_library, filename)


def process_file(gene_set_library, filename):
    gene_lists = get_gene_lists(filename)
    if not gene_lists:  # skip if the gene lists are empty
        print(f'No genes found in {filename}. Skipping...')
        return

    for gene_list_name, genes in gene_lists.items():
        if not genes:  # skip if the gene list is empty
            continue

        genes_str = '\n'.join(genes)

        # Extract the network summary type from the filename
        network_summary_type = re.search(r'(.+)_NETWORK_SUMMARY\.txt', filename).group(1)

        output_filename = os.path.join(OUTPUT_DIR, f"{network_summary_type}_{gene_list_name}_{gene_set_library}.csv")

        print(f'Processing {filename} ({gene_list_name})...')
        user_list_id = add_gene_list(genes_str)
        result = run_enrichment(user_list_id, gene_set_library)


        # Print enrichment results as a table
        print('Enrichment Results:')
        table_header = ["Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes",
                        "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
        table_data = []
        for item in result[gene_set_library]:
            rank, term_name, p_value, z_score, combined_score, overlapping_genes, adjusted_p_value, old_p_value, old_adjusted_p_value = item
            table_data.append([rank, term_name, p_value, z_score, combined_score, overlapping_genes,
                               adjusted_p_value, old_p_value, old_adjusted_p_value])

        print("\t".join(table_header))
        for i, row in enumerate(table_data):
            if i >= 5:  # Stop after printing 5 rows
                break
            print("\t".join(str(cell) for cell in row))

        print()


        # Export enrichment results to CSV
        with open(output_filename, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(table_header)
            writer.writerows(table_data)

        print(f"Enrichment results exported to {output_filename}")
        print()



def add_gene_list(genes):

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    description = 'Example gene list'

    payload = {
        'list': (None, genes),
        'description': (None, description)
    }

    while True:  # loop indefinitely until a successful request
        time.sleep(2)
        response = requests.post(ENRICHR_URL, files=payload)
        if response.status_code == 429:  # Too Many Requests
            print("Received 429 response. Sleeping and retrying.")
            time.sleep(10)  # wait 10 seconds before retrying
            continue
        elif not response.ok:  # any other unsuccessful status
            response.raise_for_status()

        data = json.loads(response.text)
        return data['userListId']


def run_enrichment(user_list_id, gene_set_library):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    while True:  # loop indefinitely until a successful request
        time.sleep(2)
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
        )
        if response.status_code == 429:  # Too Many Requests
            print("Received 429 response. Sleeping and retrying.")
            time.sleep(10)  # wait 10 seconds before retrying
            continue
        elif not response.ok:  # any other unsuccessful status
            response.raise_for_status()

        data = json.loads(response.text)
        return data


def get_gene_lists(filename):
    gene_lists = {
        'degree': [],
        'betweenness': [],
    }

    with open(filename, 'r') as file:
        content = file.readlines()

    current_section = None
    eigenvector_class = None

    for line in content:
        # Remove leading/trailing whitespace
        line = line.strip()

        # Check if we are in the Degree or Betweenness section
        if 'Node |     Degree' in line:
            current_section = 'degree'
        elif 'Node | Betweenness Centrality' in line:
            current_section = 'betweenness'
        elif 'Node | Eigenvector Centrality' in line:
            current_section = 'eigenvector'
        elif re.match(r'CLASS \d+:', line):
            # If we encounter a CLASS line in Eigenvector section, create a new list
            if current_section == 'eigenvector':
                eigenvector_class = re.search(r'CLASS (\d+):', line).group(1)
                gene_lists[f'eigenvector_class_{eigenvector_class}'] = []
        elif current_section and line and not line.startswith('-'):
            # If we're in a section, and this is a line with gene data
            # (not empty, not a '---' line), add gene to the appropriate list
            gene_match = re.match(r'\w+', line)
            if gene_match:
                gene_name = gene_match.group(0)
                if current_section != 'eigenvector':
                    gene_lists[current_section].append(gene_name)
                elif eigenvector_class is not None:
                    gene_lists[f'eigenvector_class_{eigenvector_class}'].append(gene_name)

    return gene_lists



def move_files_to_directories():
    base_directory = "./enrichment_results/"  # The directory where your files are located
    files = os.listdir(base_directory)

    # Define the subdirectories
    subdirs = ['top_degree', 'top_betweenness', 'top_eigenvector_classes']

    # Create directories if they don't exist
    for subdir in subdirs:
        if not os.path.exists(os.path.join(base_directory, subdir)):
            os.makedirs(os.path.join(base_directory, subdir))

    for file in files:
        # Only move csv files
        if file.endswith('.csv'):
            if "degree" in file:
                shutil.move(os.path.join(base_directory, file), os.path.join(base_directory, 'top_degree', file))
            elif "betweenness" in file:
                shutil.move(os.path.join(base_directory, file), os.path.join(base_directory, 'top_betweenness', file))
            elif "eigenvector" in file:
                shutil.move(os.path.join(base_directory, file),
                            os.path.join(base_directory, 'top_eigenvector_classes', file))



if __name__ == "__main__":
    gene_set_library_list = [
        "GO_Biological_Process_2021",
        "GO_Molecular_Function_2021",
        "KEGG_2021_Human",
        "PPI_Hub_Proteins",
        "DrugMatrix",
        "MAGMA_Drugs_and_Diseases",
        "IDG_Drug_Targets_2022"
    ]

    filenames = [filename for filename in os.listdir('.') if
                 re.match(r'.+(genes|interactors|omim_genes_interactions)_NETWORK_SUMMARY\.txt$',
                          filename)]

    for filename in filenames:
        for gene_set_library in gene_set_library_list:
            process_file(gene_set_library, filename)

    # Organize files into directories
    move_files_to_directories()

    print()
    print('-------------------------------------------------------------------------------------')
    print()
    print('                     Thank you for using 丂尸丫ᗪ🝗尺爪工爪    ❤')
    print()
    print('-------------------------------------------------------------------------------------')
    print()

"""

-------------------+
Gene set libraries |
-------------------+

_______________________________________________________________________________________________

Diseases / Drugs

    COVID-19_Related_Gene_Sets_2021
    Orphanet_Augmented_2021
    LINCS_L1000_Chem_Pert_Consensus_Sigs
    LINCS_L1000_CRISPR_KO_Consensus_Sigs
    GTEx_Aging_Signatures_2021
    ClinVar_2019
    HDSigDB_Human_2021
    HDSigDB_Mouse_2021
    DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019
    PheWeb_2019
    PhenGenI_Association_2021
    Proteomics_Drug_Atlas_2023
    DepMap_WG_CRISPR_Screens_Broad_CellLines_2019
    TG_GATES_2020
    GWAS_Catalog_2019
    UK_Biobank_GWAS_v1
    DisGeNET
    DSigDB
    ARCHS4_IDG_Coexp
    DrugMatrix
    Old_CMAP_up
    Old_CMAP_down
    GeneSigDB
    OMIM_Disease
    OMIM_Expanded
    VirusMINT
    MSigDB_Computational
    MSigDB_Oncogenic_Signatures
    Virus_Perturbations_from_GEO_up
    Virus_Perturbations_from_GEO_down
    dbGaP
    Tunica_Media
    Achilles_fitness_increase
    Achilles_fitness_decrease
    Rare_Diseases_AutoRIF_ARCHS4_Predictions
    Rare_Diseases_GeneRIF_ARCHS4_Predictions
    Rare_Diseases_GeneRIF_Gene_Lists
    Rare_Diseases_AutoRIF_Gene_Lists
    MAGMA_Drugs_and_Diseases
    Diabetes_Perturbations_GEO_2022
    IDG_Drug_Targets_2022

Pathways:

    Reactome_2022
    BioPlanet_2019
    WikiPathway_2021_Human
    KEGG_2021_Human
    ARCHS4_Kinases_Coexp
    Elsevier_Pathway_Collection
    MSigDB_Hallmark_2020
    BioCarta_2016
    HumanCyc_2016
    NCI-Nature_2016
    Panther_2016
    BioPlex_2017
    huMAP
    PPI_Hub_Proteins
    KEA_2015
    Kinase_Perturbations_from_GEO_down
    Kinase_Perturbations_from_GEO_up
    Virus-Host_PPI_P-HIPSTer_2020
    NURSA_Human_Endogenous_Complexome
    CORUM
    SILAC_Phosphoproteomics
    HMS_LINCS_KinomeScan
    Phosphatase_Substrates_from_DEPOD
    SubCell_BarCode
    PFOCR_Pathways
    Metabolomics_Workbench_Metabolites_2022
    GlyGen_Glycosylated_Proteins_2022
    The_Kinase_Library_2023

Ontologies:

    GO_Biological_Process_2021
    GO_Molecular_Function_2021
    GO_Cellular_Component_2021
    MGI_Mammalian_Phenotype_Level_4_2021
    Human_Phenotype_Ontology
    Jensen_TISSUES
    Jensen_COMPARTMENTS
    Jensen_DISEASES
    SynGO_2022
    KOMP2_Mouse_Phenotypes_2022

Transcription:

    ChEA_2022
    ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
    ARCHS4_TFs_Coexp
    TF_Perturbations_Followed_by_Expression
    TRRUST_Transcription_Factors_2019
    FANTOM6_lncRNA_KD_DEGs
    lncHUB_lncRNA_Co-Expression
    Enrichr_Submissions_TF-Gene_Coocurrence
    TRANSFAC_and_JASPAR_PWMs
    Epigenomics_Roadmap_HM_ChIP-seq
    TargetScan_microRNA_2017
    miRTarBase_2017
    ENCODE_TF_ChIP-seq_2015
    TF-LOF_Expression_from_GEO
    ENCODE_Histone_Modifications_2015
    Transcription_Factor_PPIs
    Genome_Browser_PWMs

"""
