import numpy as np

class Node:
    def __init__(self, name, left=None, right=None):
        self.name = name
        self.left = left 
        self.right = right
        self.parent = None
        self.transition_prob = None # transition probability matrix from node to parent
    
    def is_leaf(self):
        return(self.left is None and self.right is None)

    def is_root(self):
        return(self.parent is None)

    def compute_transition_probabilities(self, data, header):
        if not self.is_leaf():
            self.left.compute_transition_probabilities(data, header)
            self.right.compute_transition_probabilities(data, header)
        elif self.is_root():
            return
        else:
            y = data[:, header.index(self.parent.name + "freq")]
            x = data[:, header.index(self.name + "freq")]
            counts, x_, y_ = np.histogram2d(x, y, bins=np.linspace(0, 1.01, 102))
            dist = self.normalize(counts) 
            self.transition_prob = dist
            return

    def compute_likelihood(self):
        # compute P(self | parent)

        if self.left is not None:
            left_likelihood = self.left.compute_likelihood()

        if self.right is not None:
            right_likelihood = self.right.compute_likelihood()
        
        return