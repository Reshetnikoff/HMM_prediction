import pandas as pd
import numpy as np
import sys
from os.path import exists
import os
from sklearn.model_selection import StratifiedKFold
from joblib import Parallel, delayed
from config import threshold_path, output_path 

genes = os.listdir('./output_HMM') #list of .csv files with domains likelihoods
drugs = ['Isoniazid', 'Kanamycin', 'Ethambutol', 'Capreomycin', 'Ciprofloxacin', 
'Moxifloxacin', 'Ofloxacin', 'Pyrazinamide', 'Rifampicin', 'Amikacin', 'Streptomycin', 
'Prothionamide', 'Ethionamide']
first_line_drugs = ['Isoniazid', 'Ethambutol', 'Pyrazinamide', 'Rifampicin', 'Streptomycin']
phenofolder='../db/pheno' #path to files with phenotypes
source = '../db/annotated_data' #path to translations

#we get dictionary with domain likelihoods, domain name, dictionary with broken genes of each isolate, 
#list of isolates in the training sample and cooresponding list of phenotypes, and return likelihood threshold
def threshold_determination(domain_lh, domain, broken_genes, X, y):
    gene = domain[:domain.rfind('_')]
    sample_info = {}
    if 'wild' in domain_lh:
        wild_lh = float(domain_lh['wild'])
    else:
        #print('this domain has no wild type info')
        return 0
    likelihoods = []
    for sample in X:
        sample_name = sample + '_result.tsv'
        if sample_name in domain_lh:
            likelihoods.append(float(domain_lh[sample_name]))
        else:
            if gene not in broken_genes[sample_name]:
                likelihoods.append(float(wild_lh))
            else:
                likelihoods.append(0)
    likelihoods = np.array(likelihoods)
    nonbroken_ids = np.where(likelihoods != 0)
    preds, pheno = likelihoods[nonbroken_ids], y[nonbroken_ids]
    res_count = len(pheno[pheno == 1])
    sus_count = len(pheno[pheno == 0])
    if (res_count * sus_count == 0):
        return 0
    if len(preds[preds != wild_lh]) == 0:
        return 0
    res_idx, sus_idx = np.where(pheno == 1), np.where(pheno == 0)
    deltas = np.abs(preds - wild_lh)
    threshold = np.max(deltas) + 1
    r_deltas = deltas[res_idx]
    s_deltas = deltas[sus_idx]
    max_difference = 0
    for delta in np.unique(deltas):
        if delta == 0:
            continue
        fraction_s = len(s_deltas[s_deltas >= delta])/sus_count
        fraction_r = len(r_deltas[r_deltas >= delta])/res_count
        if np.abs(fraction_s - fraction_r) > max_difference:
            threshold = delta
            max_difference = np.abs(fraction_s - fraction_r)
    return threshold

#we get a list of broken genes for each isolate to exclude them when determining likelihood threshold
def get_broken_genes(folder):
    broken_genes = {}
    for f in os.listdir(folder):
        broken_genes[f] = []
        with open(folder + f) as info:
            for line in info:
                if line.rfind('broken') != -1:
                    broken_genes[f].append(line.strip().split('\t')[0])
    #print(broken_genes)
    return broken_genes

#we make startified split of sample, we write down the split to *_fold_split.txt file
def stratified_split(my_drug, outfolder):
    skf = StratifiedKFold(n_splits=5, shuffle=True) #random_state=42 if we use it only once
    X, y = [], []
    with open(phenofolder + my_drug + '.pheno') as f:
        for line in f:
            X.append(line.strip().split('\t')[0])
            y.append(int(line.strip().split('\t')[1]))
    X, y = np.array(X), np.array(y)
    train_indices = []
    c = 1
    for train_index, test_index in skf.split(X, y):
        train_indices.append(train_index)
        X_train, X_test = X[train_index], X[test_index]
        with open(outfolder + '.' + str(c) + '_fold_split.txt', 'w') as logfile:
            logfile.write('train\n')
            for sample in X_train:
                logfile.write(sample + '\n')
            logfile.write('test\n')
            for sample in X_test:
                logfile.write(sample + '\n')
        c += 1
    return X, y, train_indices

#we read files with likelihoods and make a dictionary with domain info for each isolate
def process_drug(my_drug, infolder, outfolder):
    os.makedirs(outfolder, exist_ok=True)
    broken_genes = get_broken_genes(source + my_drug + '/')
    X, y, train_indices = stratified_split(my_drug, outfolder)
    for gene in genes:
        domain_dict = {}
        with open(infolder + '/' + gene) as gene_df:
            gene_df.readline()
            for line in gene_df:
                domain, pred, sample = line.strip().split('\t')
                if domain not in domain_dict:
                    domain_dict[domain] = {}
                if sample not in domain_dict[domain]:
                    domain_dict[domain][sample] = pred
        output = open(outfolder + '/' + gene, 'w')
        output.write('name\tdrug\tt1\tt2\tt3\tt4\tt5\n')
        for domain in domain_dict:
            output.write(domain + '\t' + my_drug)
            for train_index in train_indices:
                X_train, y_train = X[train_index], y[train_index]
                t = threshold_determination(domain_dict[domain], domain, broken_genes, X_train, y_train)
                output.write('\t' + str(t))
            output.write('\n')
        output.close()
    
if __name__ == "__main__":
    infolder = './HMM_output/' #folder with domain likelihoods (should have subfolders for each drug)
    outfolder = './domains_thresholds_max_dif_2/' #folder with threshold

    try:
        os.makedirs(outfolder)
    except FileExistsError:
        pass

    os.makedirs(outfolder, exist_ok=True)
    tasks = Parallel(n_jobs=13)(delayed(process_drug)(my_drug, infolder + my_drug, 
                                                      outfolder + my_drug) for my_drug in drugs)
