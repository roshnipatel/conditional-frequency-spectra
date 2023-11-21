import sys
import fastDTWF
import numpy as np
import torch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--prefix')
parser.add_argument('--suffix')
parser.add_argument('--n_ancestral')
parser.add_argument('--h_coeff')
parser.add_argument('--s_coeff')
parser.add_argument('--out')
args = parser.parse_args()

n_ancestral = int(float(args.n_ancestral)) * 2

# Parse selection and dominance coefficients 
s = -float(args.s_coeff)
h = float(args.h_coeff)
if s == 0: # neutral
    het, hom = 0, 0
elif h == 0.5: # negative selection
    het = s * h
    hom = s
else: # stabilizing selection
    het = s * h
    hom = 0

mu = torch.tensor(1e-8, dtype=torch.float64)
zero = torch.tensor(0., dtype=torch.float64)

s_het = torch.tensor(het, dtype=torch.float64)
s_hom = torch.tensor(hom, dtype=torch.float64)

base_ps = fastDTWF.wright_fisher_ps_mutate_first(
    n_ancestral, mu, zero, s_het, s_hom
)
index_sets = fastDTWF.coarse_grain(
    base_ps, n_ancestral, 0.1, True, False
).detach().numpy()
num_sets = np.max(index_sets) + 1

to_return = None

for chunk_idx in range(21):
    start = int(num_sets / 20 * chunk_idx)
    end = int(num_sets / 20 * (chunk_idx + 1))
    end = min([end, num_sets])

    fname = args.prefix + str(chunk_idx) + args.suffix
    x = np.load(fname)
    if to_return is None:
        to_return = np.zeros((n_ancestral, x.shape[1]))
    for i in range(start, end):
        if index_sets[-1] == i:
            continue
        keep = (index_sets == i)[:-1]
        to_return[keep] = x[i-start]

np.save(args.out, to_return)
