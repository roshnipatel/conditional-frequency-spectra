"""
Author: Jeffrey P. Spence
"""
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--nested')
parser.add_argument('--disjoint')
parser.add_argument('--out')
args = parser.parse_args() 

nested_sim = np.load(args.nested)
disjoint_sim = np.load(args.disjoint)

ld_nested_sim = np.zeros((199, 199))
ld_disjoint_sim = np.zeros((199, 199))

for i in range(199):
    for j in range(199):
        x = np.zeros(200)
        y = np.zeros(200)
        x[0:(i+1)] = 1
        y[0:(j+1)] = 1
        ld_nested_sim[i, j] = np.corrcoef(x, y)[0, 1]

        if (i+1) + (j+1) <= 200:
            x = np.zeros(200)
            y = np.zeros(200)
            x[0:(i+1)] = 1
            y[(i+1):(i+1+j+1)] = 1
            ld_disjoint_sim[i, j] = np.corrcoef(x, y)[0, 1]

ld_sq_nested_sim = ld_nested_sim**2
ld_sq_disjoint_sim = ld_disjoint_sim**2

p_derived_linked = np.zeros((11, 199))
for t_idx, t in enumerate([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999]):
    print(nested_sim[ld_sq_nested_sim >= t].sum() / (nested_sim[ld_sq_nested_sim >= t].sum() + disjoint_sim[ld_sq_disjoint_sim >= t].sum()))
    for i in range(199):
        keep_nested = ld_sq_nested_sim[i] >= t
        keep_disjoint = ld_sq_disjoint_sim[i] >= t
        p_derived_linked[t_idx, i] = nested_sim[i, keep_nested].sum() / (nested_sim[i, keep_nested].sum() + disjoint_sim[i, keep_disjoint].sum())

np.savetxt(args.out + "prob_derived_linked.txt", p_derived_linked)

freqs_1 = []
freqs_2 = []
for i in range(199):
    for k, count in zip(np.where(ld_sq_nested_sim[i] >= 0.5)[0], nested_sim[i, ld_sq_nested_sim[i] >= 0.5]):
        count = int(count)
        freqs_2.extend([(k+1)/200] * count)
        freqs_1.extend([(i+1)/200]*count)
    for k, count in zip(np.where(ld_sq_disjoint_sim[i] >= 0.5)[0], disjoint_sim[i, ld_sq_disjoint_sim[i] >= 0.5]):
        count = int(count)
        freqs_2.extend([1.-(k+1)/200] * count)
        freqs_1.extend([(i+1)/200]*count)
freqs_1 = np.array(freqs_1)
freqs_2 = np.array(freqs_2)
freqs_1[freqs_1 > 0.5] = 1. - freqs_1[freqs_1 > 0.5]
freqs_2[freqs_2 > 0.5] = 1. - freqs_2[freqs_2 > 0.5]

np.savetxt(args.out + "linked_freqs_1.txt", freqs_1)
np.savetxt(args.out + "linked_freqs_2.txt", freqs_2)
