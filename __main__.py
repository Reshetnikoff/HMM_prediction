import sys
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
from config import NUM_THREADS

if __name__ == "__main__":
    start_time = time.time()

    # Get input_file
    samples_file = sys.argv[1]
    samples_file = f"./sample_batches/{samples_file}"
    output = f"./HMM_output/{samples_file}/"

    samples_list = pd.read_csv(
        samples_file,

        header=None,
        index_col=None,
        sep='\t',
        squeeze=True,
        names=["samples"]
    )
    samples_list = samples_list.to_numpy()

    # Get support files block
    f = open(gene_file, 'r')
    genes = [x[:-1] for x in f.readlines()][:-1]
    f.close()

    HMMmodel = collections.namedtuple('HMMmodel', ["name", "model"])
    annotation = get_annotation(file_annotation)
    reference = get_reference(file_reference)
    pfam = get_pfam(pfam_file)

    domains = np.unique([x.split("_")[0] for x in pfam['domain']])

    # Get HMM models block
    HMMmodels = []
    print("Create HMM profiles")
    get_HMM_profile(domains)
    pool = Pool(NUM_THREADS)
    results = pool.map(get_HMM_profile, get_batches(domains, NUM_THREADS))
    for value in results:
        HMMmodels.extend(value)
    pool.close()
    pool.join()
    print("HMM profiles have created!")
    print("--- %s seconds ---" % (time.time() - start_time))

    print(len(HMMmodels))

    # Get data block
    drug = collections.namedtuple('drug', ['samples', 'data'])
    data = []
    for sample_file in samples_list:
        temp_data = pd.read_csv(

            f"{data_path}{sample_file}",

            sep='\t', header=None, index_col=None,
            names=['gene', 'pos', "ind1", "ind2", "act"],
            keep_default_na=False, na_values=['']
        )
        data.append(temp_data)
    samples_data = drug(samples=samples_list[0].values, data=data)

    # Main block for prediction
    procs = []
    for batch in get_batches(genes, NUM_THREADS):
        proc = Process(target=get_score, args=(batch, pfam, annotation, reference,
                                               samples_data, HMMmodels, "./output"))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()
    print(f"List of samples: {samples_list} is processed")
    print("--- %s seconds ---" % (time.time() - start_time))