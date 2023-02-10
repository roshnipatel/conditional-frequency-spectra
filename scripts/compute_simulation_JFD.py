import numpy as np
from jfd import SimulatedDist
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--simulated_data', nargs='+')
parser.add_argument('--ancestral_dist')
parser.add_argument('--out')
args = parser.parse_args()

simulated = SimulatedDist()
for fp in args.simulated_data:
    simulated.add_observations(filepath=fp)
simulated.initialize()

ancestral_dist = np.loadtxt(args.ancestral_dist, skiprows = 1, usecols = 1)
counts, _ = np.histogram(ancestral_dist, bins=simulated.freq_bin_bounds)
ancestral_dist = simulated.normalize(counts)
simulated.set_ancestral_distribution(ancestral_dist)

simulated.compute_ancestral_conditional_distribution("CHBfreq")
simulated.compute_ancestral_conditional_distribution("YRIfreq")

simulated.compute_marginal_study_distribution()

simulated.compute_study_conditional_distribution("ancestral_freq")

YRI = simulated.compute_study_conditional_distribution("YRIfreq")
np.savetxt(args.out + "YRI.txt", YRI, delimiter='\t')
CHB = simulated.compute_study_conditional_distribution("CHBfreq")
np.savetxt(args.out + "CHB.txt", CHB, delimiter='\t')

simulated.ancestral_pop = "CEUfreq" # weird hack to compute the dist I want in the next line
dist = simulated.compute_ancestral_conditional_distribution("ooAfreq")
np.savetxt(args.out + "ooA.txt", dist, delimiter='\t')