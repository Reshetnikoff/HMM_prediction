# PFAM domains likelihood prediction and strong changes in domains detection

This repository contains code to obtain PFAM domains likelihoods based on lists of mutations observed in Mycobacterium tuberculosis (MTB) isolate genome. In addition, domains likelihoods can be easily binarized to point out protein domains that have been altered the most.

## Usage

First of all, you should make configuration changes in the config.py file. You should specify the following paths:

 - **HMM_file**, a path to the PFAM HMM file
 - **pfam_file**, a path to domains info file. Columns from left to right represent domain name, domain end position within a gene, corresponding gene name, and domain start position within a gene
 - **data_path**, a path to directory with MTB isolates mutation data (detailed information on files format is given below)
 - **output_path**, a path to directory for domains likelihoods predictions
 - **threshold_path**, a path to directory with likelihood thresholds that will be used for binarization
 - **bin_path**, a path to binarized likelihoods, will contain a list of domains for each MTB isolate that were strongly altered by mutations
 - **file_annotation**, a path to annotation of M.Tuberculosis genome file
 - **file_reference**, a path to MTB H37Rv reference file  
 - **gene_file**, list of MTB genes considered
 - **NUM_THREADS**, a number of threads

Isolate mutations data file must contain five tab-separated columns. First column has a gene name. Second column has position of a mutation. Next columns depend on type of mutations. If it's SNP third, fourth and fifth contain wild type aminoacide or nucleotide, mutated aminoacide or nucleotide and 'snp' symbol. If it's aminoacid insertion the columns contains wild type aminoacids in the position, inserted sequence and 'ins' symbol. If it's nucleotide insertion third column has dash unlike the aninoacid insertion. If it's animoacid deletion the columns have deleted sequence, wild type sequence in the position and 'del' sign. Nucleotide deletion columns have only a number of deleted position and 'del' symbol. If it's left clip mutation only third column represent 'left clip symbol'. If it's right clip third and fourth columns has a number of clipped position and 'right_clip' symbol. Third and fourth right extension columns have extented sequence and 'right_extension symbol'. Third and fourth left extension columns have extented sequence and 'left_extension symbol'. Also our test dataset has 'broken gene' features but they are ignored.

In addition, you should download PFAM HMM file:

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz 

First, the code allows to get a list of likelihoods predictions for all available Pfam domains and for every isolate. You can run program by command:

        python make_prediction.py

Output directory will contain a file with predictions for every gene for every isolate in the input directory. Each file will contain three columns. The columns represent from left to right domain name, logarithmic prediction and name of file or 'wild' symbol.

Next script defines a likelihood threshold for each drug that is further can be used for binarization:

	python domains_threshold_determination_folds.py

And finally domain_feature_generation_folds.py can be used to obtain binarized data:
	
	python domain_feature_generation_folds.py

