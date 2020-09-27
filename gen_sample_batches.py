import os

import math

N = 5

input_dir = "../new_data"
output_dir = "./sample_batches"

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

samples_list = os.listdir(input_dir)
samples_list.remove("all_sample_data.csv")

for ind, samples in enumerate(get_batches(samples_list, 5)):
    f = open(f"{output_dir}/{ind}.csv", "w")
    samples = "\n".join(samples)
    f.write(samples)
    f.close()
  

