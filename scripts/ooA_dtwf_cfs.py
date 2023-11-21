import numpy as np
from node import Node
from helpers import compute_joint_probability, bin_ancestral_dist, bin_dtwf_based_transition
from probability_dist import ProbabilityDist, Axis
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--positive', action="store_true")
parser.add_argument('--ancestral')
parser.add_argument('--ancestor_yri')
parser.add_argument('--ancestor_ooa')
parser.add_argument('--ooa_ceu')
parser.add_argument('--ooa_chb')
parser.add_argument('--out_prefix')
args = parser.parse_args()

ceu = Node("CEUfreq")
chb = Node("CHBfreq")
yri = Node("YRIfreq")
ooa = Node("ooAfreq", left=ceu, right=chb)
ancestor = Node("ancestral_freq", left=yri, right=ooa)

ancestral_sfs = np.load(args.ancestral)
binned_ancestral_sfs, ancestral_bins = bin_ancestral_dist(ancestral_sfs, condition_segregating=False)
ancestor.marginal = ProbabilityDist(binned_ancestral_sfs, [Axis(ancestor.name, ancestral_bins)])

# forward matrices
yri_forward = np.load(args.ancestor_yri).T
yri_forward = np.append(yri_forward, np.fmax(np.array([1 - np.sum(yri_forward, axis=0)]), 0), axis=0)
if args.positive:
    yri_forward = np.rot90(np.rot90(yri_forward))
binned_dist, yri_bins, _ = bin_dtwf_based_transition(yri_forward)
yri_conditional_ancestor = ProbabilityDist(binned_dist, [Axis(yri.name, yri_bins), Axis(ancestor.name, ancestral_bins)],
                                          conditional_var=Axis(ancestor.name, ancestral_bins))
yri.transition = yri_conditional_ancestor
yri_conditional_ancestor.write(args.out_prefix + "yri_conditional_ancestor.txt")

ooa_forward = np.load(args.ancestor_ooa).T
ooa_forward = np.append(ooa_forward, np.fmax(np.array([1 - np.sum(ooa_forward, axis=0)]), 0), axis=0)
if args.positive:
    ooa_forward = np.rot90(np.rot90(ooa_forward))
binned_dist, ooa_bins, _ = bin_dtwf_based_transition(ooa_forward)
ooa_conditional_ancestor = ProbabilityDist(binned_dist, [Axis(ooa.name, ooa_bins), Axis(ancestor.name, ancestral_bins)],
                                          conditional_var=Axis(ancestor.name, ancestral_bins))
ooa.transition = ooa_conditional_ancestor
ooa_conditional_ancestor.write(args.out_prefix + "ooa_conditional_ancestor.txt")

ceu_forward = np.load(args.ooa_ceu).T
ceu_forward = np.append(ceu_forward, np.fmax(np.array([1 - np.sum(ceu_forward, axis=0)]), 0), axis=0)
if args.positive:
    ceu_forward = np.rot90(np.rot90(ceu_forward))
binned_dist, ceu_bins, _ = bin_dtwf_based_transition(ceu_forward)
ceu_conditional_ooa = ProbabilityDist(binned_dist, [Axis(ceu.name, ceu_bins), Axis(ooa.name, ooa_bins)],
                                          conditional_var=Axis(ooa.name, ooa_bins))
ceu.transition = ceu_conditional_ooa
ceu_conditional_ooa.write(args.out_prefix + "ceu_conditional_ooa.txt")

chb_forward = np.load(args.ooa_chb).T
chb_forward = np.append(chb_forward, np.fmax(np.array([1 - np.sum(chb_forward, axis=0)]), 0), axis=0)
if args.positive:
    chb_forward = np.rot90(np.rot90(chb_forward))
binned_dist, chb_bins, _ = bin_dtwf_based_transition(chb_forward)
chb_conditional_ooa = ProbabilityDist(binned_dist, [Axis(chb.name, chb_bins), Axis(ooa.name, ooa_bins)],
                                          conditional_var=Axis(ooa.name, ooa_bins))
chb.transition = chb_conditional_ooa
chb_conditional_ooa.write(args.out_prefix + "chb_conditional_ooa.txt")

# joint matrices
joint_ancestor_children, children_conditional_ancestor = compute_joint_probability(ancestor) # P(ancestor, children)
joint_children = joint_ancestor_children.marginalize(ancestor.name) # P(children)

# remaining forward matrices
drop_yri = children_conditional_ancestor.marginalize(yri.name) # P(ceu, chb | ancestor)

ceu_conditional_ancestor = drop_yri.marginalize(chb.name) # P(ceu | ancestor)
ceu_conditional_ancestor.write(args.out_prefix + "ceu_conditional_ancestor.txt")

chb_conditional_ancestor = drop_yri.marginalize(ceu.name) # P(chb | ancestor)
chb_conditional_ancestor.write(args.out_prefix + "chb_conditional_ancestor.txt")

# backward matrices
drop_chb = joint_ancestor_children.marginalize(chb.name)

joint_ancestor_ceu = drop_chb.marginalize(yri.name) # P(ancestor, ceu)
ancestor_conditional_ceu, ceu_marg = joint_ancestor_ceu.condition(ceu.name) # P(ancestor | ceu)
ancestor_conditional_ceu.write(args.out_prefix + "ancestor_conditional_ceu.txt")

ooa_marginal = ooa.transition.multiply_by_marginal(ancestor.marginal).marginalize(ancestor.name)
joint_ooa_ceu = ceu.transition.multiply_by_marginal(ooa_marginal)
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