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
from config import data_path, output_path
from config import NUM_THREADS

if __name__ == "__main__":
    start_time = time.time()

    drug = sys.argv[1]
    data_path += f'/{drug}'
    output_path += f'/{drug}'

    # Get input_file
    samples_list = os.listdir(data_path)
    samples_list.remove(drug)

    # Config output
    try:
        os.makedirs(output_path)
    except FileExistsError:
        pass

    # Get support files block
    f = open(gene_file, 'r')
    genes = [x[:-1] for x in f.readlines()][:-1]
    f.close()

    HMMmodel = collections.namedtuple('HMMmodel', ["name", "model"])
    annotation = get_annotation(file_annotation)
    reference = get_reference(file_reference)
    pfam = get_pfam(pfam_file)

    domains = ['PF00204', 'PF00986', 'PF01751', 'PF02518']

    # Get HMM models block
    HMMmodels = []
    print("Create HMM profiles")
    HMMmodels = get_HMM_profile(domains)
    print("HMM profiles have created!")
    print("--- %s seconds ---" % (time.time() - start_time))


    # Get data block
    drug = collections.namedtuple('drug', ['samples', 'data'])
    data = []
    for sample_file in samples_list:
        temp_data = pd.read_csv(

            f"{data_path}/{sample_file}",

            sep='\t', header=None, index_col=None,
            names=['gene', 'pos', "ind1", "ind2", "act"],
            keep_default_na=False, na_values=[''], usecols=[0, 1, 2, 3, 4]
        )
        data.append(temp_data)
    samples_data = drug(samples=samples_list, data=data)

    # Main block for prediction

    gene = 'gyrB'
    get_score([gene], pfam, annotation, reference,
                samples_data, HMMmodels, output_path)
    print(f"List of samples: {data_path} is processed and write in {output_path}")
    print("--- %s seconds ---" % (time.time() - start_time))