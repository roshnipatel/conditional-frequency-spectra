import numpy as np
from node import Node
from helpers import compute_sim_based_transition, compute_joint_probability, bin_ancestral_dist
from probability_dist import ProbabilityDist, Axis
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--simulated_data", nargs='+')
parser.add_argument('--ancestral_sfs')
parser.add_argument('--out_prefix')
args = parser.parse_args()

ancestral_bins = [0, 0.00003, 0.0001, 0.0005, 0.001, 0.002, 
                  0.005, 0.01, 0.015, 0.02, 0.025] + \
                 [i / 100 for i in range(3, 100)] + [1]
ancestral_sfs = np.load(args.ancestral_sfs)
binned_ancestral_sfs, _ = bin_ancestral_dist(ancestral_sfs, bins=ancestral_bins, condition_segregating=True)

ceu = Node("CEUfreq", pop_size=31721)
chb = Node("CHBfreq", pop_size=65653)
yri = Node("YRIfreq", pop_size=3721)
ooa = Node("ooAfreq", pop_size=2831, left=ceu, right=chb)
ancestor = Node("ancestral_freq", pop_size=23721, left=yri, right=ooa)
ancestor.marginal = ProbabilityDist(binned_ancestral_sfs, [Axis(ancestor.name, ancestral_bins)])

# load simulations from file to generate transition matrices for tree
sims = None
for fp in args.simulated_data:
    with open(fp, 'r') as f:
        header = f.readline().strip().replace('"', '').split('\t')
    new_sims = np.loadtxt(fp, skiprows = 1)
    if sims is None:
        sims = new_sims
    else:
        sims = np.concatenate([sims, new_sims], axis=0)

compute_sim_based_transition(ancestor, sims, header, ancestral_bins)

joint_ancestor_children, children_conditional_ancestor = compute_joint_probability(ancestor) # P(ancestor, children)
joint_children = joint_ancestor_children.marginalize(ancestor.name) # P(children)

# forward matrices
yri_conditional_ancestor = yri.transition # P(yri | ancestor)
yri_conditional_ancestor.write(args.out_prefix + "yri_conditional_ancestor.txt")

drop_yri = children_conditional_ancestor.marginalize(yri.name) # P(ceu, chb | ancestor)

ceu_conditional_ancestor = drop_yri.marginalize(chb.name) # P(ceu | ancestor)
ceu_conditional_ancestor.write(args.out_prefix + "ceu_conditional_ancestor.txt")

ooa_marginal = ooa.transition.multiply_by_marginal(ancestor.marginal).marginalize(ancestor.name)
joint_ooa_ceu = ceu.transition.multiply_by_marginal(ooa_marginal)
ceu_conditional_ooa, _ = joint_ooa_ceu.condition(ooa.name)
ceu_conditional_ooa.write(args.out_prefix + "ceu_conditional_ooa.txt")

chb_conditional_ancestor = drop_yri.marginalize(ceu.name) # P(chb | ancestor)
chb_conditional_ancestor.write(args.out_prefix + "chb_conditional_ancestor.txt")

# backward matrices
drop_chb = joint_ancestor_children.marginalize(chb.name)

joint_ancestor_ceu = drop_chb.marginalize(yri.name) # P(ancestor, ceu)
ancestor_conditional_ceu, ceu_marg = joint_ancestor_ceu.condition(ceu.name) # P(ancestor | ceu)
ancestor_conditional_ceu.write(args.out_prefix + "ancestor_conditional_ceu.txt")

ooa_conditional_ceu, _ = joint_ooa_ceu.condition(ceu.name)
ooa_conditional_ceu.write(args.out_prefix + "ooa_conditional_ceu.txt")

joint_ancestor_yri = drop_chb.marginalize(ceu.name) # P(ancestor, yri)
ancestor_conditional_yri, yri_marg = joint_ancestor_yri.condition(yri.name) # P(ancestor | yri)
ancestor_conditional_yri.write(args.out_prefix + "ancestor_conditional_yri.txt")

# marginal matrices
ooa_marginal.write(args.out_prefix + "ooa_marginal.txt")
ceu_marg.write(args.out_prefix + "ceu_marginal.txt")
yri_marg.write(args.out_prefix + "yri_marginal.txt")

# conditional matrices
modern_conditional_ceu, _ = joint_children.condition(ceu.name) # P (chb, yri | ceu)
modern_conditional_ceu.write(args.out_prefix + "modern_conditional_ceu.txt")

yri_ceu = joint_children.marginalize(chb.name)
yri_chb = joint_children.marginalize(ceu.name)
chb_ceu = joint_children.marginalize(yri.name)

yri_conditional_ceu, _ = yri_ceu.condition(ceu.name) # P(yri | ceu)
yri_conditional_ceu.write(args.out_prefix + "yri_conditional_ceu.txt")

chb_conditional_ceu, _ = chb_ceu.condition(ceu.name) # P(chb | ceu)
chb_conditional_ceu.write(args.out_prefix + "chb_conditional_ceu.txt")

ceu_conditional_yri, _ = yri_ceu.condition(yri.name) # P(ceu | yri)
ceu_conditional_yri.write(args.out_prefix + "ceu_conditional_yri.txt")

chb_conditional_yri, _ = yri_chb.condition(yri.name) # P(chb | yri)
chb_conditional_yri.write(args.out_prefix + "chb_conditional_yri.txt")

ceu_conditional_chb, _ = chb_ceu.condition(chb.name) # P(ceu | chb)
ceu_conditional_chb.write(args.out_prefix + "ceu_conditional_chb.txt")

yri_conditional_chb, _ = yri_chb.condition(chb.name) # P(yri | chb)
yri_conditional_chb.write(args.out_prefix + "yri_conditional_chb.txt")