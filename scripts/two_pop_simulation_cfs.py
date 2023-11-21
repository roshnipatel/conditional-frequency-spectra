import numpy as np
from node import Node
from helpers import compute_sim_based_transition, compute_joint_probability, bin_ancestral_dist
from probability_dist import ProbabilityDist, Axis
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--p1_sims", nargs='+')
parser.add_argument("--p2_sims", nargs='+')
parser.add_argument("--n_p1")
parser.add_argument("--n_p2")
parser.add_argument('--ancestral_sfs')
parser.add_argument('--out_prefix')
args = parser.parse_args()

ancestral_bins = [0, 0.0001, 0.0005, 0.001, 0.002,
                  0.005, 0.01, 0.015, 0.02, 0.025] + \
                 [i / 100 for i in range(3, 100)] + [1]
ancestral_sfs = np.load(args.ancestral_sfs)
binned_ancestral_sfs, _ = bin_ancestral_dist(ancestral_sfs, bins=ancestral_bins, condition_segregating=True)

p1 = Node("p1freq", pop_size=int(float(args.n_p1)))
p2 = Node("p2freq", pop_size=int(float(args.n_p2)))
ancestor = Node("ancestral_freq", pop_size=ancestral_sfs.shape[0], left=p1, right=p2)
ancestor.marginal = ProbabilityDist(binned_ancestral_sfs, [Axis(ancestor.name, ancestral_bins)])

# load simulations from file to generate transition matrices for tree
p1_sims = None
with open(args.p1_sims[0], 'r') as f:
    p1_header = f.readline().strip().replace('"', '').split('\t')
for fp in args.p1_sims:
    new_p1_sims = np.loadtxt(fp, skiprows = 1)
    if p1_sims is None:
        p1_sims = new_p1_sims
    else:
        p1_sims = np.concatenate([p1_sims, new_p1_sims], axis=0)
p1_header[p1_header.index("modern_freq")] = p1.name
compute_sim_based_transition(p1, p1_sims, p1_header, ancestral_bins)

p2_sims = None
with open(args.p2_sims[0], 'r') as f:
    p2_header = f.readline().strip().replace('"', '').split('\t')
for fp in args.p2_sims:
    new_p2_sims = np.loadtxt(fp, skiprows = 1)
    if p2_sims is None:
        p2_sims = new_p2_sims
    else:
        p2_sims = np.concatenate([p2_sims, new_p2_sims], axis=0)
p2_header[p2_header.index("modern_freq")] = p2.name
compute_sim_based_transition(p2, p2_sims, p2_header, ancestral_bins)

joint_ancestor_children, children_conditional_ancestor = compute_joint_probability(ancestor) # P(ancestor, children)
joint_children = joint_ancestor_children.marginalize(ancestor.name) # P(children)

# forward matrices
p1_conditional_ancestor = p1.transition # P(p1 | ancestor)
p1_conditional_ancestor.write(args.out_prefix + "p1_conditional_ancestor.txt")

p2_conditional_ancestor = p2.transition # P(p2 | ancestor)
p2_conditional_ancestor.write(args.out_prefix + "p2_conditional_ancestor.txt")

# backward matrices
joint_ancestor_p1 = joint_ancestor_children.marginalize(p2.name) # P(ancestor, p1)
ancestor_conditional_p1, p1_marg = joint_ancestor_p1.condition(p1.name) # P(ancestor | p1)
ancestor_conditional_p1.write(args.out_prefix + "ancestor_conditional_p1.txt")
p1_marg.write(args.out_prefix + "p1_marginal.txt")

joint_ancestor_p2 = joint_ancestor_children.marginalize(p1.name) # P(ancestor, p2)
ancestor_conditional_p2, p2_marg = joint_ancestor_p2.condition(p2.name) # P(ancestor | p2)
ancestor_conditional_p2.write(args.out_prefix + "ancestor_conditional_p2.txt")
p2_marg.write(args.out_prefix + "p2_marginal.txt")

# conditional matrices
p1_conditional_p2, _ = joint_children.condition(p2.name) # P(p1 | p2)
p1_conditional_p2.write(args.out_prefix + "p1_conditional_p2.txt")

p2_conditional_p1, _ = joint_children.condition(p1.name) # P(p2 | p1)
p2_conditional_p1.write(args.out_prefix + "p2_conditional_p1.txt")