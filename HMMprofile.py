import numpy as np
import pandas as pd
import os
import time
import math
from functools import partial


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable

import multiprocessing
from multiprocessing import Process
import collections

from pomegranate import *


HMM_file = "../db/Pfam-A.hmm"
pfam_file = "../db/PfamDomainsInfo.csv"
data_path = "../new_data/"
pheno = "../new_pheno/"
file_annotation = "../db/AL123456_rev.gff"
file_reference = '../db/h37rv.fasta'
GENE_FILE = "./genes.csv"
NUM_THREADS = 30



class profileHMM():
    def __init__(self, domain, HMM_file, verbose=False):
        self.HMM_file = HMM_file
        self.domain = domain
        self.size = 0
        self.model = None
        self.error = 0
        
        
        self.m_trans_prob = []
        self.i_trans_prob = []
        self.d_trans_prob = []
        self.St_initial_prob = []
        self.I0_initial_prob = []
        self.D_initial_prob = []
        
        data = self.parse_HMM(HMM_file)
        if len(data) == 0:
            self.error = 1
        else:
            self.amino_acids = data[0].split()[1:]

            trans_prob, self.emission_prob_match, self.emission_prob_ins = self.get_prob(data)
            self.transform_trans_prob(trans_prob)

            if verbose == True:
                print(len(self.emission_prob_match))
                print(len(self.emission_prob_ins))
                print(size)

    def predict(self, seq):
        return(self.model.log_probability(list(seq)))
    
    def transform_trans_prob(self, trans_prob):
        self.m_trans_prob.append([[x[0], x[1], x[2]] for x in trans_prob[1:]])
        self.i_trans_prob.append([[x[3], x[4]] for x in trans_prob[1:]])
        self.d_trans_prob.append([[x[5], x[6]] for x in trans_prob[1:]])
        self.St_initial_prob.append(trans_prob[0][:3])
        self.I0_initial_prob.append(trans_prob[0][3:5])

        self.m_trans_prob = np.array(self.m_trans_prob[0], dtype=np.float64)
        self.i_trans_prob = np.array(self.i_trans_prob[0], dtype=np.float64)
        self.d_trans_prob = np.array(self.d_trans_prob[0], dtype=np.float64)
        self.St_initial_prob = np.array(self.St_initial_prob, dtype=np.float64)
        self.I0_initial_prob = np.array(self.I0_initial_prob, dtype=np.float64)



        self.d_trans_prob = self.round_num(self.d_trans_prob)
        self.i_trans_prob = self.round_num(self.i_trans_prob)
        self.m_trans_prob = self.round_num(self.m_trans_prob)
        self.St_initial_prob = self.round_num(self.St_initial_prob)
        self.I0_initial_prob = self.round_num(self.I0_initial_prob)
        
    def get_prob(self, data):
        trans_prob = []
        emission_prob_match = []
        emission_prob_ins = []
        
        for ind in range(int((len(data)-2)/3)):
            trans_prob.append(data[(ind+1)*3+1].split()) 
            emission_prob_ins.append(data[(ind+1)*3].split()) 
            if ind == 0:
                continue
            else:
                emission_prob_match.append(data[(ind+1)*3-1].split()[1:-5])


        trans_prob = np.array(trans_prob)
        emission_prob_match = np.array(emission_prob_match)
        emission_prob_ins = np.array(emission_prob_ins)

        trans_prob[trans_prob == "*"] = np.inf
        emission_prob_match[emission_prob_match == "*"] = np.inf
        emission_prob_ins[emission_prob_ins == "*"] = np.inf

        trans_prob = np.array(trans_prob, dtype=np.float64)
        emission_prob_match = np.array(emission_prob_match, dtype=np.float64)
        emission_prob_ins = np.array(emission_prob_ins, dtype=np.float64)

        emission_prob_match = self.round_num(emission_prob_match)
        emission_prob_ins = self.round_num(emission_prob_ins)
        
        return(trans_prob, emission_prob_match, emission_prob_ins)
        
    
    def round_num(self, arr):
        new_arr = []
        for x in arr:
            new_x = np.around(np.exp(-x), decimals=5)
            diff = (1.0 - np.sum(new_x))
            new_arr.append([new_x[0] + diff, *new_x[1:]])
        new_arr = np.array(new_arr, dtype=np.float64)
        new_arr = np.around(new_arr, decimals=5)
        return(new_arr)

    def parse_HMM(self, HMM_file):
        data_file = open(self.HMM_file, 'r')

        data = []

        in_hmm = False
        save = False
        for line in data_file.readlines():
            if self.domain in line:
                in_hmm = True
            if in_hmm == True:
                if line[:4] == "LENG":
                    self.size = int(line.split(" ")[-1])
                if line[:3] == "HMM":
                    save = True
            if save == True:
                if line[:2] == "//":
                    break
                data.append(line)

        data_file.close() 
        return data
    
    def build(self):
        model = HiddenMarkovModel()

        I0 = State(DiscreteDistribution({x : y for (x, y) in zip(self.amino_acids, self.emission_prob_ins[0])}))
        M1 = State(DiscreteDistribution({x : y for (x, y) in zip(self.amino_acids, self.emission_prob_match[0])}))
        D1 = State(None)

        model.add_states(I0, M1, D1)

        model.add_transition(model.start, M1, self.St_initial_prob[0][1])
        model.add_transition(model.start, I0, self.St_initial_prob[0][0])
        model.add_transition(model.start, D1, self.St_initial_prob[0][2])

        model.add_transition(I0, M1, self.I0_initial_prob[0][0])
        model.add_transition(I0, I0, self.I0_initial_prob[0][1])
        


        for i in range(1, self.size):
            if i == 1:
                Ii = I0
                Mi = M1
                Di = D1
            else:
                Ii = Iii
                Mi = Mii
                Di = Dii

            Iii = State(DiscreteDistribution({x : y for (x, y) in zip(self.amino_acids, self.emission_prob_ins[i])}))
            Mii = State(DiscreteDistribution({x : y for (x, y) in zip(self.amino_acids, self.emission_prob_match[i])}))
            Dii = State(None)
            
            model.add_states(Iii, Mii, Dii)

            model.add_transition(Mi, Mii, self.m_trans_prob[i-1][0])
            model.add_transition(Mi, Iii, self.m_trans_prob[i-1][1])
            model.add_transition(Mi, Dii, self.m_trans_prob[i-1][2])

            model.add_transition(Iii, Mii, self.i_trans_prob[i-1][0])
            model.add_transition(Iii, Iii, self.i_trans_prob[i-1][1])

            model.add_transition(Di, Mii, self.d_trans_prob[i-1][0])
            model.add_transition(Di, Dii, self.d_trans_prob[i-1][1])

        Ii = Iii
        Iii = State(DiscreteDistribution({x : y for (x, y) in zip(self.amino_acids, self.emission_prob_ins[-1])}))
        
        model.add_states(Iii)
        
        model.add_transition(Mii, model.end, self.m_trans_prob[-1][0])
        model.add_transition(Mii, Iii, self.m_trans_prob[-1][1])

        model.add_transition(Iii, model.end, self.i_trans_prob[-1][0])
        model.add_transition(Iii, Iii, self.i_trans_prob[-1][1])

        model.bake()
        
        self.model = model


def round_num(arr):
    new_arr = []
    for x in arr:
        new_x = np.around(np.exp(-x), decimals=5)
        diff = (1.0 - np.sum(new_x))
        new_arr.append([new_x[0] + diff, *new_x[1:]])
    new_arr = np.array(new_arr, dtype=np.float64)
    new_arr = np.around(new_arr, decimals=5)
    return(new_arr)


def change_seq(seq, data, diff=0):
    shift = 0
    len_seq = len(seq)
    for (gene, pos, ind1, ind2, act) in data.values:
        pos = pos - diff + 1
        if act == "snp":
            if pos+shift >= len(seq):
                seq = seq[:pos-1+shift] + ind2
            else:
                if ind2 == "-":
                    seq = seq[:pos-1+shift] + seq[pos+shift:]
                else:
                    seq = seq[:pos-1+shift] + ind2 + seq[pos+shift:]
        if act == "ins":
            # интересная вещь с этой вставкой - ограничиться длиной или нет? - открытый вопрос
            # такого if по идее не должно произойти
            if pos+shift-1 >= len(seq):
                continue
            seq = seq[:pos-1+shift] + ind2 + seq[pos-1+shift:]
            shift += len(ind2)
        if act == "del":
            if pos-1+len(ind1)+shift >= len(seq):
                continue
            seq = seq[:pos-1+shift] + seq[pos-1+len(ind1)+shift:]
            shift -= len(ind1)
        if ind1 == "left_clip":
            shift -= pos
            seq = seq[pos:]
        if ind2 == "right_clip":
            ind1 = int(ind1)
            seq = seq[:pos-1]
        if ind2 == "right_extension":
            seq = seq[:pos-1] + ind1[:len_seq-pos+1]
        if ind2 == "left_extension":
            seq = ind1[-pos:] + seq[pos:]
    return seq

def delete_last_space(s):
    if s[-1] == " ":
        return s[:-1]
    else:
        return s

def get_annotation(file):
    annotation = pd.read_csv(file, sep='\t', header=None, index_col=None)
    annotation = annotation[annotation[2] == 'Gene']
    annotation.drop(columns=[0, 1, 2, 5, 7], inplace=True)

    annotation[10] = ([x.split(" ")[-1] for x in annotation[8]])

    temp = ([delete_last_space(x.split(";")[0]) for x in annotation[8]])
    annotation[5] = ([str(x).split(" ")[-1] for x in temp])

    bi_genes = pd.DataFrame()
    bi_genes["gene"] = annotation[10]
    bi_genes["alt_gene"] = annotation[5]
    bi_genes.to_csv("bi_genes.csv", sep='\t', header=True, index=None)

    annotation.drop(columns=[10 ,8], inplace=True)
    annotation.columns = ['start', 'end', 'strand','gene']
    return annotation


def get_reference(file):
    reference_file = open(file)
    reference = [x[:-1] for x in reference_file.readlines()[1:]] 
    reference = "".join(reference)
    reference_file.close()
    return reference


def get_gene_seq(gene, annotation, reference):
    st, end, strand, _ = annotation[annotation["gene"] == gene].values[0]
    seq = Seq(reference[st-1:end], generic_dna)
    if strand == '-':
        seq = seq.reverse_complement()
    elif strand == '.' and seq.__str__()[:3] != "ATG":
        seq = seq.reverse_complement()
    seq = seq.translate(table=CodonTable.unambiguous_dna_by_id[3])
    seq = seq.__str__()[:-1]
    return(seq)


def get_pfam(pfam_file):
    return(pd.read_csv(pfam_file, sep='\t', index_col=False, header=0))


def check_model(HMMmodels, domain):
    for name, model in HMMmodels:
        if domain == name:
            return model
    return -1

def get_score(genes):
    for gene in genes:
        preds = pd.DataFrame(columns=["name", "pred", "sample"])    
        domains = pfam[pfam['gene'] == gene]

        seq = get_gene_seq(gene, annotation, reference)
        for (domain, end, gene, start) in domains.values:
            clear_domain = domain.split("_")[0]
            result = check_model(HMMmodels, clear_domain)
            if result == -1:
                continue
            model = result
            
            domain_seq = seq[start-1:end]
            wild_pred = model.predict(domain_seq)
            preds = preds.append({"name" : f"{gene}_{domain}", "pred" : wild_pred, "sample" : "wild"}, ignore_index=True)



        for (sample, sample_data) in zip(drugs.samples, drugs.data):
            gene_data = sample_data[sample_data['gene'] == gene]
            gene_data = gene_data[gene_data["pos"] != "changed"]
            gene_data = gene_data[gene_data["pos"] != "broken"]
            gene_data = gene_data.astype({"pos" : np.int32})
            gene_data.sort_values(by="pos", axis=0, inplace=True)

            for (domain, end, gene, start) in domains.values:
                domain_data = gene_data[gene_data['pos'] >= int(start)]
                domain_data = domain_data[gene_data['pos'] <= int(end)]
                if len(domain_data) == 0:
                    continue
                clear_domain = domain.split("_")[0]
            
                # Не очень красивый код для ускорения
                result = check_model(HMMmodels, clear_domain)
                if result == -1:
                    continue
                model = result
              
                domain_seq = seq[start-1:end]
                # Проверить change_seq
                mut_seq = change_seq(domain_seq, domain_data, diff=start)
                mut_pred = model.predict(mut_seq)
                preds = preds.append({"name" : f"{gene}_{domain}", "pred" : mut_pred, "sample" : sample}, ignore_index=True)
        preds.to_csv(f"{output}{gene}.csv", sep='\t', index=False, header=True)

def get_HMM_profile(domains):
    HMMmodels_func = []
    for domain in domains:
        model = profileHMM(domain, HMM_file)
        if model.error != 1:
             model.build()
             HMMmodels_func.append(HMMmodel(name=domain, model=model))
    return HMMmodels_func   

def get_batches(arr, n):
    batches = []
    size = len(arr)
    batch_size = math.floor(size / n)
    rest = len(arr) % n
    st = 0
    end = batch_size
    for i in range(n):
        if rest != 0:
            end += 1
            rest -= 1
        batches.append(list(arr[st:end]))
        st = end
        end += batch_size
    return batches

if __name__ == "__main__":
    # Костыль
    samples_ind = sys.argv[1]
    samples_list = f"./sample_batches/{samples_ind}.csv"
    print(samples_list)
    output = f"./HMM_output/{samples_ind}/"


    start_time = time.time()
    f = open(GENE_FILE, 'r')
    genes = [x[:-1] for x in f.readlines()][:-1]    
    f.close()
   
    HMMmodel = collections.namedtuple('HMMmodel', ["name", "model"])    
    annotation = get_annotation(file_annotation)
    reference = get_reference(file_reference)
    pfam = get_pfam(pfam_file)

    domains = np.unique([x.split("_")[0] for x in pfam['domain']])

    HMMmodels = []
    print("Create HMM profiles")
    get_HMM_profile(domains)
    # pool = multiprocessing.Pool(NUM_THREADS)
    # results = pool.map(get_HMM_profile, get_batches(domains, NUM_THREADS))
    # for value in results:
    #     HMMmodels.extend(value)
    # pool.close()
    # pool.join()
    print("HMM profiles have created!")
    print("--- %s seconds ---" % (time.time() - start_time))

    print(len(HMMmodels))
    drug = collections.namedtuple('drug', ['samples', 'data'])
    drug_data = pd.read_csv(
                                f"{samples_list}", 

                                header=None, 
                                index_col=None, 
                                sep='\t'
                            )

    data = []
    for sample_file in drug_data[0].values:
        temp_data = pd.read_csv(

               f"{data_path}{sample_file}",

               sep='\t', header=None, index_col=None,
               names=['gene', 'pos', "ind1", "ind2", "act"],
               keep_default_na=False, na_values=['']
               )
        data.append(temp_data)


    drugs = drug(samples=drug_data[0].values, data=data)
    procs = []
    for batch in get_batches(genes, NUM_THREADS):  
        proc = Process(target=get_score, args=(batch, ))   
        procs.append(proc)
        proc.start()
    
    for proc in procs:
        proc.join()
    print(f"List of samples: {samples_list} is processed")
    print("--- %s seconds ---" % (time.time() - start_time))



