import numpy as np
from node import Node
from helpers import compute_joint_probability, bin_ancestral_dist, bin_dtwf_based_transition
from probability_dist import ProbabilityDist, Axis
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--positive', action="store_true")
parser.add_argument('--ancestral')
parser.add_argument('--p1_forward')
parser.add_argument('--p2_forward')
parser.add_argument('--out_prefix')
args = parser.parse_args()

p1 = Node("p1freq")
p2 = Node("p2freq")
ancestor = Node("ancestral_freq", left=p1, right=p2)

ancestral_sfs = np.load(args.ancestral)
binned_ancestral_sfs, ancestral_bins = bin_ancestral_dist(ancestral_sfs, condition_segregating=False)
ancestor.marginal = ProbabilityDist(binned_ancestral_sfs, [Axis(ancestor.name, ancestral_bins)])

# forward matrices
p1_forward = np.load(args.p1_forward).T
p1_forward = np.append(p1_forward, np.fmax(np.array([1 - np.sum(p1_forward, axis=0)]), 0), axis=0)
if args.positive:
    p1_forward = np.rot90(np.rot90(p1_forward))
binned_dist, p1_bins, _ = bin_dtwf_based_transition(p1_forward)
p1_conditional_ancestor = ProbabilityDist(binned_dist, [Axis(p1.name, p1_bins), Axis(ancestor.name, ancestral_bins)],
                                          conditional_var=Axis(ancestor.name, ancestral_bins))
p1.transition = p1_conditional_ancestor
p1_conditional_ancestor.write(args.out_prefix + "p1_conditional_ancestor.txt")

p2_forward = np.load(args.p2_forward).T
p2_forward = np.append(p2_forward, np.fmax(np.array([1 - np.sum(p2_forward, axis=0)]), 0), axis=0)
if args.positive:
    p2_forward = np.rot90(np.rot90(p2_forward))
binned_dist, p2_bins, _ = bin_dtwf_based_transition(p2_forward)
p2_conditional_ancestor = ProbabilityDist(binned_dist, [Axis(p2.name, p2_bins), Axis(ancestor.name, ancestral_bins)],
                                          conditional_var=Axis(ancestor.name, ancestral_bins))
p2.transition = p2_conditional_ancestor
p2_conditional_ancestor.write(args.out_prefix + "p2_conditional_ancestor.txt")

# joint matrices
joint_ancestor_children, children_conditional_ancestor = compute_joint_probability(ancestor) # P(ancestor, children)
joint_children = joint_ancestor_children.marginalize(ancestor.name) # P(children)

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