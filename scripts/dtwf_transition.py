import fastDTWF.fastDTWF as fastDTWF
import numpy as np
import torch
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ancestral')
parser.add_argument('--n_ancestral')
parser.add_argument('--n_modern')
parser.add_argument('--n_gen')
parser.add_argument('--chunk')
parser.add_argument('--h_coeff')
parser.add_argument('--s_coeff')
parser.add_argument('--out')
args = parser.parse_args()

# Specify ancestral Ne
n_haploid_ancestral = int(float(args.n_ancestral)) * 2

# Parse modern Ne
n_haploid_modern = int(float(args.n_modern)) * 2

# Parse generations 
num_gen = int(float(args.n_gen))

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

chunk_idx = int(args.chunk)
assert chunk_idx >= 0
assert chunk_idx <= 20

# Load ancestral SFS
anc = np.zeros(n_haploid_ancestral + 1)
anc[:-1] = np.load(args.ancestral)
anc = torch.tensor(anc, dtype=torch.float64)

def get_likelihoods(N, s_het, s_hom):
    anc_p = anc

    s_het = torch.tensor(s_het, dtype=torch.float64)
    s_hom = torch.tensor(s_hom, dtype=torch.float64)

    base_ps = fastDTWF.wright_fisher_ps_mutate_first(
        n_haploid_ancestral, mu, zero, s_het, s_hom
    )
    index_sets = fastDTWF.coarse_grain(
        base_ps, n_haploid_ancestral, 0.1, True, False
    ).detach().numpy()
    num_sets = np.max(index_sets) + 1

    start = int(num_sets / 20 * chunk_idx)
    end = int(num_sets / 20 * (chunk_idx + 1))
    end = min([end, num_sets])

    to_return = np.zeros((end-start, N))

    for i in range(start, end):
        if i == index_sets[-1]:
            continue
        p_vec = torch.zeros(n_haploid_ancestral + 1, dtype=torch.float64)
        p_vec[index_sets == i] = anc_p[index_sets == i]
        if p_vec.sum() <= 0:
            p_vec[index_sets == i] = 1.
        p_vec /= p_vec.sum()
        p_vec = fastDTWF._integrate_likelihood_constant_size(
            vec=p_vec,
            interval_pop_size=N,
            interval_length=5,
            s_het=s_het,
            mu_0_to_1=mu,
            mu_1_to_0=zero,
            dtwf_tv_sd=0.1,
            dtwf_row_eps=1e-8,
            no_fix=False,
            sfs=False,
            injection_rate=0.,
            s_hom=s_hom,
            use_condensed=False
        )
        p_vec = fastDTWF._integrate_likelihood_constant_size(
            vec=p_vec,
            interval_pop_size=N,
            interval_length=95,
            s_het=s_het,
            mu_0_to_1=mu,
            mu_1_to_0=zero,
            dtwf_tv_sd=0.1,
            dtwf_row_eps=1e-8,
            no_fix=False,
            sfs=False,
            injection_rate=0.,
            s_hom=s_hom,
            use_condensed=True,
            refresh_gens=25
        )
        p_vec = fastDTWF._integrate_likelihood_constant_size(
            vec=p_vec,
            interval_pop_size=N,
            interval_length=num_gen-100,
            s_het=s_het,
            mu_0_to_1=mu,
            mu_1_to_0=zero,
            dtwf_tv_sd=0.1,
            dtwf_row_eps=1e-8,
            no_fix=False,
            sfs=False,
            injection_rate=0.,
            s_hom=s_hom,
            use_condensed=True,
            refresh_gens=200
        )
        to_return[i-start] = p_vec.detach().numpy()[:-1]

    return to_return

xy = get_likelihoods(n_haploid_modern, het, hom)
np.save(args.out, xy)
