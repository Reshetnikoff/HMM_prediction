# HMM prediction

The programm is used to get preidiction by PFAM's hidden markov model for domains based on list of mutations. 

## Usage

First of all, it's needed to make own configuration of the config.py file:

 - HMM_file - is a path to the PFAM HMM file
 - pfam_file - is a path to domains info file. Columns from left to right represen domain name, end position of the domain within a gene, gene name, and start position of the domain within the gene
 - data_path - is path to a dir which contain sample's data
 - output_path - is path to a output dir
 - file_annotation - is a path to annotation of M.Tuberculosis genome file
 - file_reference - is a path to reference file  
 - gene_file - is a path to list of genes
 - NUM_THREADS - a number of threads

Input file must constain five columns. First column has gene name. Second column has position of a mutation. Next columns depend on type of mutations. If it's SNP third, fourth and fifth contain wild type aminoacide or nucleotide, mutated aminoacide or nucleotide and 'snp' symbol. If it's aminoacid insertion the columns contains wild type aminoacids in the position, inserted sequence and 'ins' symbol. If it's nucleotide insertion third column has dash unlike the aninoacid insertion. If it's animoacid deletion the columns have deleted sequence, wild type sequence in the position and 'del' sign. Nucleotide deletion columns have only a number of deleted position and 'del' symbol. If it's left clip mutation only third column represent 'left clip symbol'. If it's right clip third and fourth columns has a number of clipped position and 'right_clip' symbol. Third and fourth right extension columns have extented sequence and 'right_extension symbol'. Third and fourth left extension columns have extented sequence and 'left_extension symbol'. Also our test dataset has 'broken gene' features but they are ignored.

Output directory constains list of file with prediction for every gene for every samples in the input directory. Each file constain three columns. The columns represent from left to right domain name, logarithmic prediction and name of file or 'wild' symbol.

Run program by command:

	python __main__.py  
