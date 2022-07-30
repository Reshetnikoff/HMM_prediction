import pandas as pd
import numpy as np
import sys
from os.path import exists
import os
from joblib import Parallel, delayed
from config import threshold_path, bin_path

#input_folder = '/export/data/kchukreev/domains_thresholds_max_dif_2/'
#output_folder = '/export/data/kchukreev/domains_features_max_dif_2/'
drugs = ['Isoniazid', 'Kanamycin', 'Ethambutol', 'Capreomycin', 
'Ciprofloxacin', 'Moxifloxacin', 'Ofloxacin', 'Pyrazinamide', 
'Rifampicin', 'Amikacin', 'Streptomycin', 'Prothionamide', 
'Ethionamide']
first_line_drugs = ['Isoniazid', 'Ethambutol', 'Pyrazinamide', 'Rifampicin', 'Streptomycin']
#source = '/export/data/kchukreev/HMM_031121/translation_leaves_190221_filtered_constant_upstream_intergenic_v2/' #path to folder with domains likelihoods
source = '/export/data/kchukreev/HMM_23072022/translation_leaves_190221_filtered_constant_upstream_intergenic_fixed/'
phenofolder = '/export/data/kchukreev/pheno_combined/' #path to files with phenotypes
genes = os.listdir('./output_test') #list of .csv files with domains likelihoods

#generate domain features
def generate_feature(my_drug, additional_samples, input_folder, output_folder):
    end_of_s = '_result.tsv'
    phenotypes = {line.strip().split('\t')[0] + end_of_s: 
                  int(line.strip().split('\t')[1]) for line in open(phenofolder + my_drug + '.pheno').readlines()}
    samples_domains = {}
    for gene in genes:
        domain_thresholds = {}
        domain_dict = {}
        with open(input_folder + '/' + gene) as f:
            f.readline()
            for line in f:
                temp = line.strip().split('\t')
                domain, drug = temp[0], temp[1]
                thresholds = [float(t) for t in temp[2:]]
                domain_thresholds[domain] = thresholds
        if additional_samples:
            gene_source = source + gene
        else:
            gene_source = source + my_drug + '/' + gene
        with open(gene_source) as gene_df:
            gene_df.readline()
            for line in gene_df:
                domain, pred, sample = line.strip().split('\t')
                if domain not in domain_thresholds:
                    continue
                if domain not in domain_dict:
                    domain_dict[domain] = {}
                if sample not in domain_dict[domain]:
                    domain_dict[domain][sample] = float(pred)
        for domain in domain_dict:
            if 'wild' in domain_dict[domain]:
                wild_lh = float(domain_dict[domain]['wild'])
                for sample in domain_dict[domain]:
                    if sample not in samples_domains:
                        samples_domains[sample] = []
                    sample_delta = np.abs(domain_dict[domain][sample] - wild_lh)
                    c = 0
                    for threshold in domain_thresholds[domain]:
                        c += 1
                        if (threshold != 0) and (sample_delta >= threshold):
                            samples_domains[sample].append((domain, c))
            else:
                continue
    print("I made samples_domains dict! Let me write it down")
    for sample in samples_domains:
        if (sample in phenotypes) or ((my_drug not in first_line_drugs) and (sample[:-len(end_of_s)] in additional_samples)):
            folds_changes = {}
            for change in samples_domains[sample]:
                domain, fold = change
                if fold not in folds_changes:
                    folds_changes[fold] = []
                folds_changes[fold].append(domain)
            for fold in folds_changes:
                with open(output_folder + '/fold' + str(fold) + '/' + sample, 'a') as output:
                    for domain in folds_changes[fold]:
                        output.write(domain + '\tchanged\n')

if __name__ == "__main__":
    input_folder = threshold_path #folder with thresholds
    output_folder = bin_path

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    with_additional_samples = False #True if we consider samples susceptible to first-line drugs to be also susceptible to second-line drugs
    additional_samples = []
    if with_additional_samples:
        source = source[:-1] + '_merged/'
        first_line_susceptible = {}
        for drug in first_line_drugs:
            if drug != 'Streptomycin':
                first_line_susceptible[drug] = []
                with open(phenofolder + drug + '.pheno') as f:
                    for line in f:
                        if line.strip().split('\t')[1] == '0':
                            first_line_susceptible[drug].append(line.strip().split('\t')[0])
        susc_lists = [set(first_line_susceptible[drug]) for drug in first_line_susceptible]
        additional_samples = set.intersection(*susc_lists)
    for drug in drugs:
        os.makedirs(output_folder + drug, exist_ok=True)
        for i in range(1, 6):
            os.makedirs(output_folder + drug + '/fold' + str(i), exist_ok=True)
    tasks = Parallel(n_jobs=14)(delayed(generate_feature)(my_drug, additional_samples, input_folder + my_drug, output_folder + my_drug) 
                                                          for my_drug in drugs)
