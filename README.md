# PFAM domains likelihood prediction and strong changes in domains detection

This repository contains code to obtain PFAM domains likelihoods based on lists of mutations observed in  *Mycobacterium tuberculosis* (MTB) isolates genomes. In addition, domains likelihoods can be easily binarized to point out protein domains that have been altered the most.

## Usage

First of all, you should make configuration changes in the config.py file. You should specify the following paths:

 - **HMM_file**, a path to the PFAM HMM file
 - **pfam_file**, a path to domains info file. Columns from left to right represent domain name, domain end position within a gene, corresponding gene name, and domain start position within a gene
 - **data_path**, a path to directory with MTB isolates mutation data (detailed information on files format is given below)
 - **output_path**, a path to directory for domains likelihoods predictions (an output of make_prediction.py will be written there)
 - **threshold_path**, a path to directory for likelihood thresholds that will be used for binarization (an output of domains_threshold_determination_folds.py will be written here)
 - **bin_path**, a path to binarized likelihoods, will contain a list of domains for each MTB isolate that are strongly altered by mutations (an output of domain_feature_generation_folds.py will be written there)
 - **file_annotation**, a path to annotation of M.Tuberculosis genome file
 - **file_reference**, a path to MTB H37Rv reference file  
 - **gene_file**, list of MTB genes considered
 - **NUM_THREADS**, a number of threads

Isolate mutations data file has mutation annotation. Files of a proper format can be obtained with Annotation.py script from the following repository: https://github.com/dashabykova/MTB_project. Annotation file contains five tab-separated columns: gene name, a mutation position according to corresponding protein sequence, reference amino acid, alternative amino acid and mutation type. Mutation type can be either 'snp' for single amino acid substitution, 'ins' for an insertion and 'del' for a deletion. For insertions there is a dash ('-') symbol for reference amino acid sequence, for deletions there is no reference or alternative amino acids, just a deletion length (in amino acids) instead. Also, there can be additional labels in annotation such as 'right_extension' or 'left_extension' if a protein sequence is prolonged, and 'right_clip' and 'left_clip' if a protein sequence is shortened due to start or stop codon mutation. If a protein sequence is shortened by more than a half or prolonged by more than 30%, the gene is considered to be 'broken' and it is ignored during domains likelihood calculation.         

In addition, you should download PFAM HMM file:

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz 

The following command allows to get a list of likelihoods predictions for all available PFAM domains for every isolate in a data_path folder:

        python make_prediction.py

An output directory will contain a file with predictions for every gene from a gene_file. Each output file contains information on all gene's PFAM domains alterations in all isolates from a data_path folder. The columns represent domain name, logarithmic prediction and isolate name if a domain is altered in this isolate or 'wild' label for reference protein domain likelihood.

Next script defines a likelihood threshold for each drug that can be further used for binarization (see paper for more details on a threshold determination):

	python domains_threshold_determination_folds.py


And finally domain_feature_generation_folds.py can be used to obtain binarized data:
	
	python domain_feature_generation_folds.py

An output directory will contain a file for each isolate with a list of domains that are strongly altered in this isolate. In output files each row will contain a domain name and 'changed' label.
