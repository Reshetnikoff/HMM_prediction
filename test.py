import sys
import os
import time
from multiprocessing import Process, Pool

import pandas as pd
import numpy as np

import collections

from utils import get_score
from utils import get_annotation, get_reference, get_pfam, get_HMM_profile
from utils import get_batches

from config import file_annotation, file_reference, pfam_file, gene_file
from config import data_path


if __name__ == "__main__":
    start_time = time.time()

    # Get support files block
    samples_list = os.listdir(data_path)
    f = open(gene_file, 'r')
    genes = [x[:-1] for x in f.readlines()][:-1]
    f.close()

    HMMmodel = collections.namedtuple('HMMmodel', ["name", "model"])
    annotation = get_annotation(file_annotation)
    reference = get_reference(file_reference)
    pfam = get_pfam(pfam_file)

    # Get HMMmodels block
    domain = ["PF00141"]

    HMMmodels = []
    print("Create HMM profiles")

    HMMmodels.extend(get_HMM_profile(domain))
    print(len(HMMmodels))
    print("HMM have been created")

    # Get data block
    print("Load sample's data...")
    drug = collections.namedtuple('drug', ['samples', 'data'])

    data = []
    for sample_file in samples_list:
        temp_data = pd.read_csv(

            f"{data_path}/{sample_file}",

            sep='\t', header=None, index_col=None,
            names=['gene', 'pos', "ind1", "ind2", "act"],
            keep_default_na=False, na_values=['']
        )
        data.append(temp_data)
    drug_data = drug(samples=samples_list, data=data)
    print("Data have been loaded")

    # Process data
    print("Process data...")
    get_score(["katG"], pfam, annotation, reference, drug_data, HMMmodels, "./output")
    print("Data have been processed")

    print("--- %s seconds ---" % (time.time() - start_time))