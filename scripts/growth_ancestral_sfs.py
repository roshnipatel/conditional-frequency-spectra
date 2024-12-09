import fastDTWF.fastDTWF as fastDTWF
import numpy as np
import torch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--h_coeff')
parser.add_argument('--s_coeff')
parser.add_argument('--out')
args = parser.parse_args()

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
    pop_size_list=[20000, 55000, 150000],
    switch_points=[2000, 1000, 0],
    sample_size=20000,
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
