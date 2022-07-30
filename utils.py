import math

import pandas as pd
import numpy as np
import collections

from Bio.Seq import Seq
from Bio.Data import CodonTable

from HMMprofile import profileHMM

from config import HMM_file

"""
Support HMM block
"""

def get_score(genes, pfam, annotation, reference, drugs, HMMmodels, output):
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

            domain_seq = seq[start - 1:end]
            wild_pred = model.predict(domain_seq)
            preds = preds.append({"name": f"{gene}_{domain}", "pred": wild_pred, "sample": "wild"}, ignore_index=True)


        for (sample, sample_data) in zip(drugs.samples, drugs.data):
            gene_data = sample_data[sample_data['gene'] == gene]
            gene_data = gene_data[gene_data["pos"] != "changed"]
            gene_data = gene_data[gene_data["pos"] != "broken"]
            gene_data = gene_data.astype({"pos": np.int32})
            gene_data.sort_values(by="pos", axis=0, inplace=True)


            for (domain, end, gene, start) in domains.values:
                left_gene_data = gene_data[gene_data['ind1'] == 'left_clip']
                right_gene_data = gene_data[gene_data['ind2'] == 'right_clip']

                left_gene_data = left_gene_data[left_gene_data['pos'] >= end]
                right_gene_data = right_gene_data[right_gene_data['pos'] <= start]
                if len(left_gene_data) > 0 or len(right_gene_data) > 0:
                    preds = preds.append({"name": f"{gene}_{domain}", "pred": -np.inf, "sample": sample},
                                         ignore_index=True)
                    continue


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

                domain_seq = seq[start - 1:end]
                # Проверить change_seq
                mut_seq = change_seq(domain_seq, domain_data, diff=start)
                mut_pred = model.predict(mut_seq)
                preds = preds.append({"name": f"{gene}_{domain}", "pred": mut_pred, "sample": sample},
                                     ignore_index=True)
        preds.to_csv(f"{output}/{gene}.csv", sep='\t', index=False, header=True)


def check_model(HMMmodels, domain):
    for name, model in HMMmodels:
        if domain == name:
            return model
    return -1

def get_HMM_profile(domains):
    HMMmodel = collections.namedtuple('HMMmodel', ["name", "model"])

    HMMmodels_func = []
    for domain in domains:
        model = profileHMM(domain, HMM_file)
        if model.error != 1:
             model.build()
             HMMmodels_func.append(HMMmodel(name=domain, model=model))
    return HMMmodels_func

"""
    Get function
"""

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
            seq = seq[:pos+shift] + seq[pos-1+len(ind1)+shift:]
            shift -= len(ind1) + 1
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

def get_annotation(file):
    annotation = pd.read_csv(file, sep='\t', header=None, index_col=None)
    annotation = annotation[annotation[2] == 'CDS']
    annotation.drop(columns=[0, 1, 2, 5, 7], inplace=True)

    annotation[10] = ([x.split(" ")[-1] for x in annotation[8]])

    temp = ([delete_last_space(x.split(";")[2]) for x in annotation[8]])
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
    seq = Seq(reference[st-1:end])
    if strand == '-':
        seq = seq.reverse_complement()
    elif strand == '.' and seq.__str__()[:3] != "ATG":
        seq = seq.reverse_complement()
    seq = seq.translate(table=CodonTable.unambiguous_dna_by_id[3])
    seq = seq.__str__()[:-1]
    return(seq)


def get_pfam(pfam_file):
    return(pd.read_csv(pfam_file, sep='\t', index_col=False, header=0))


"""
    Support function
"""

def delete_last_space(s):
    if s[-1] == " ":
        return s[:-1]
    else:
        return s

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

