import numpy as np
from probability_dist import ProbabilityDist, Axis

def create_bins(n_alleles):
    n_intermediate_bins = int((n_alleles - 25) / 
                              (.01 * n_alleles))
    bin_list = [0, 1, 2, 5, 10, 20] + \
        [20 + i * .01 * n_alleles for i in range(1, n_intermediate_bins + 1)] + \
        [n_alleles - 5, n_alleles - 2, n_alleles - 1, n_alleles] 
    bin_list = [i / n_alleles for i in bin_list]
    bin_list[1] = bin_list[1] + 0.000001 # heuristic to make sure allele counts of 1 end up in their own bin
    return(bin_list)            

def compute_sim_based_transition(node, data, header, ancestral_bins):
    if node.is_parent():
        if node.left is not None:
            compute_sim_based_transition(node.left, data, header, ancestral_bins)

        if node.right is not None:
            compute_sim_based_transition(node.right, data, header, ancestral_bins)

    if not node.is_ancestor():
        x = data[:, header.index(node.name)]
        y = data[:, header.index(node.parent.name)]

        # computes bin IDs for each frequency. bin IDs range from 0 (corresponding 
        # to [0]) to (n_bins - 1) (corresponding to (0.99999, 1]). functionally, the first 
        # and last bins should only contain alleles that are lost and reach fixation
        # respectively.
        x_bins = create_bins(node.pop_size * 2)
        if node.parent.is_ancestor():
            y_bins = ancestral_bins
        else:
            y_bins = create_bins(node.parent.pop_size * 2)
        n_x_bins, n_y_bins = len(x_bins), len(y_bins)
        binned_x = np.digitize(x, x_bins, right=True)
        binned_y = np.digitize(y, y_bins, right=True)

        assert np.all(binned_x < n_x_bins)
        assert np.all(binned_y < n_y_bins)
        assert np.all(binned_x >= 0)
        assert np.all(binned_y >= 0)

        # generates a histogram based on the bin ID computed by np.digitize in the
        # previous step. np.histogram2d bins these bin IDs into [0, 1), [1, 2), ...
        # [n_bins - 1, n_bins] such that again, the first and last bins only
        # contain alleles that are lost and reach fixation respectively.
        counts, x_, y_ = np.histogram2d(binned_x, binned_y, 
                                        bins=[np.array(range(n_x_bins + 1)), np.array(range(n_y_bins + 1))])

        # adjust transition matrices such that P(x_modern = 1 | x_ancestor = 1) = 1.
        # (this will not be automatically inferred from simulations.)
        counts[:, -1] = 0
        counts[-1, -1] = 1

        dist = counts / np.sum(counts, axis=0)
        node.transition = ProbabilityDist(dist, [Axis(node.name, x_bins), Axis(node.parent.name, y_bins)],
                                          conditional_var=Axis(node.parent.name, y_bins))
        return

def compute_joint_probability(node):
    if not node.is_parent():
        return(node.transition)
    
    # Compute P(children | self)
    if node.left is not None and node.right is not None:
        left = compute_joint_probability(node.left)
        right = compute_joint_probability(node.right)
        node_conditional = left.joint(right)
    elif node.left is not None:
        node_conditional = compute_joint_probability(node.left)
    elif node.right is not None:
        node_conditional = compute_joint_probability(node.right)
    else:
        raise Exception("help! node is a parent but has no children")

    if not node.is_ancestor():
        node_joint = node_conditional.multiply_by_marginal(node.transition) # P(children, self | parent)
        node_marginalized = node_joint.marginalize(node.name) # P(children | parent)
        return(node_marginalized)
    else:
        if node.marginal is None:
            raise Exception("cannot compute joint distribution without marginal distribution for ancestor")
        
        node_joint = node_conditional.multiply_by_marginal(node.marginal) # P(children, ancestor)
        return(node_joint, node_conditional)
        
def bin_ancestral_dist(count_probs, bins=None, condition_segregating=False):
    if condition_segregating:
        # condition on alleles segregating in the ancestor
        count_probs[0] = 0
        count_probs /= np.sum(count_probs)

    n_sites = count_probs.shape[0]
    if not bins:
        bins = create_bins(n_sites)
    n_bins = len(bins)
    mask = np.digitize(np.array(range(n_sites)) / n_sites, bins, right=True)
    freq_probs = np.zeros(n_bins)
    for i in range(n_bins):
        freq_probs[i] = np.sum(count_probs[mask == i])

    return(freq_probs, bins)

def bin_dtwf_based_transition(arr, condition_segregating=False):
    """Convert count-based probability matrices into corresponding frequency bins."""
    # note that columns (ancestral population) exclude fixation while rows (modern 
    # population) include fixation
    n_row, n_col = arr.shape
    row_bins = create_bins(n_row - 1)
    col_bins = create_bins(n_col)
    n_row_bins, n_col_bins = len(row_bins), len(col_bins)
    row_mask = np.digitize(np.array(range(n_row)) / (n_row - 1), row_bins, right=True)
    col_mask = np.digitize(np.array(range(n_col)) / (n_col), col_bins, right=True)
    
    tmp = np.zeros((n_row_bins, n_col))
    for i in range(n_row_bins):
        tmp[i] = np.sum(arr[row_mask == i,:], axis=0)
    
    final = np.zeros((n_col_bins, n_row_bins))
    for i in range(n_col_bins):
        final[i] = np.sum(tmp[:,col_mask == i], axis=1)
    final = final.T

    # adjust transition matrices such that P(x_modern = 1 | x_ancestor = 1) = 1
    final[n_row_bins - 1, n_col_bins - 1] = 1

    final /= np.sum(final, axis=0) 

    if condition_segregating:
        final[:,0] = 0
    
    return(final, row_bins, col_bins) 