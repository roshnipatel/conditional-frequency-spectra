import moments
from probability_dist import ProbabilityDist, Axis
import numpy as np
import argparse

# not set up to run with snakemake!
# run this in data/distributions/moments-jouganous_w_migration/

parser = argparse.ArgumentParser()
parser.add_argument('--s')
parser.add_argument('--h')
args = parser.parse_args()

# unscaled haploid population sizes
n_YRI = 23721 * 2
n_ooA = 2831 * 2
n_CEU = 2512 * 2
n_CHB = 1019 * 2

# unscaled migration rates
m_YRI_ooA = 0.000158
m_YRI_CEU = 0.000011
m_YRI_CHB = 0.0000048
m_CEU_CHB = 0.0000419

# unscaled generation times
t_epoch_1 = 2517 # (119k years - 46k years) / 29 (years/gen)
t_epoch_2 = 1586 # 46k years / 29 (years / gen)

# unscaled final haploid population sizes
n_CEU_final = n_CEU * (1 + 0.0016) ** t_epoch_2
n_CHB_final = n_CHB * (1 + 0.0026) ** t_epoch_2

# specify Ne
Ne = n_YRI

# sample size (arbitrarily taking 100 from each population)
n_samp = 100

# convert selection coefficient
gamma = 2 * Ne * float(args.s)
h = float(args.h)

# initialize frequency spectrum in ancestral population
fs = moments.LinearSystem_1D.steady_state_1D(n_samp * 3, gamma=gamma, h=h)
fs = moments.Spectrum(fs)

# split out-of-Africa branch from YRI
fs = fs.split(0, n_samp, n_samp * 2)

# scale for first epoch
n_YRI_scaled = n_YRI / Ne
n_ooA_scaled = n_ooA / Ne
t_epoch_1_scaled = t_epoch_1 / (2 * Ne)
M_epoch_1 = np.array([
    [0, m_YRI_ooA],
    [m_YRI_ooA, 0]
])
M_epoch_1_scaled = M_epoch_1 * 2 * Ne

# integrate first epoch
fs.integrate([n_YRI_scaled, n_ooA_scaled], t_epoch_1_scaled, 
             m=M_epoch_1_scaled, gamma=gamma, h=h)

# split CEU and CHB
fs = fs.split(1, n_samp, n_samp)

# scale for second epoch
n_CEU_scaled = n_CEU / Ne
n_CHB_scaled = n_CHB / Ne
n_CEU_final_scaled = n_CEU_final / Ne
n_CHB_final_scaled = n_CHB_final / Ne
t_epoch_2_scaled = t_epoch_2 / (2 * Ne)
M_epoch_2 = np.array([
    [0, m_YRI_CEU, m_YRI_CHB],
    [m_YRI_CEU, 0, m_CEU_CHB],
    [m_YRI_CHB, m_CEU_CHB, 0]
])
M_epoch_2_scaled = M_epoch_2 * 2 * Ne

# define growth function for second epoch
growth_func = lambda t: [
        n_YRI_scaled,
        n_CEU_scaled * np.exp(np.log(n_CEU_final_scaled / n_CEU_scaled) * t / t_epoch_2_scaled),
        n_CHB_scaled * np.exp(np.log(n_CHB_final_scaled / n_CHB_scaled) * t / t_epoch_2_scaled)
    ]

# integrate second epoch
fs.integrate(growth_func, t_epoch_2_scaled, 
             m=M_epoch_2_scaled, gamma=gamma, h=h)

# convert to a probability distribution
fs_scaled = fs.data / np.sum(fs.data) # moments masks non-segregating sites but we need to keep these
bins = [i / 100 for i in range(0, 101)]
joint = ProbabilityDist(fs_scaled, [Axis("YRIfreq", bins), Axis("CEUfreq", bins), Axis("CHBfreq", bins)])
ceu_conditional, _ = joint.condition("CEUfreq")
chb_ceu = ceu_conditional.marginalize("YRIfreq")
chb_ceu.write("h{0}_s{1}_chb_conditional_ceu.txt".format(args.h, args.s))
yri_ceu = ceu_conditional.marginalize("CHBfreq")
yri_ceu.write("h{0}_s{1}_yri_conditional_ceu.txt".format(args.h, args.s))