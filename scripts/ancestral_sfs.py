import fastDTWF.fastDTWF as fastDTWF
import numpy as np
import torch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sim')
parser.add_argument('--h_coeff')
parser.add_argument('--s_coeff')
parser.add_argument('--out')
args = parser.parse_args()

# Parse simulation scheme to determine ancestral Ne
if args.sim[:9] == "jouganous":
    n_diploid = 23721
elif args.sim == "gutenkunst_wo_migration":
    n_diploid = 12300
elif args.sim == "ragsdale_w_migration":
    n_diploid = 14780
else:
    n_diploid = int(float(args.sim[9:]))
n_haploid = n_diploid * 2

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

likelihoods = fastDTWF.get_likelihood(
    pop_size_list=[n_haploid],
    switch_points=[0],
    sample_size=n_haploid,
    s_het=s_het,
    s_hom=s_hom,
    mu_0_to_1=mu,
    mu_1_to_0=zero,
    dtwf_tv_sd=0.1,
    dtwf_row_eps=1e-8,
    sampling_tv_sd=0.05,
    sampling_row_eps=1e-8,
    no_fix=True,
    sfs=False,
    high_precision_stationary=True,
    injection_rate=0.
)

count_probs = likelihoods.detach().numpy()[:-1]
np.save(args.out, count_probs)
