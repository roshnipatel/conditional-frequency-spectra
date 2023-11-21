import numpy as np

class Axis:
    def __init__(self, name, bins):
        self.name = name

        assert isinstance(bins, list)
        for i in bins:
            assert isinstance(i, (int, float))

        self.bins = bins
        self.n_bins = len(bins)

    def __eq__(self, other):
        if isinstance(other, Axis):
            return((self.name == other.name) and (self.bins == other.bins))
        return NotImplemented

    def __hash__(self):
        return(hash((self.name, self.n_bins)))

    def __str__(self):
        return(self.name + ": " + ",".join([str(i) for i in self.bins])) 

class ProbabilityDist:
    def __init__(self, dist, axes, conditional_var=None):
        dist = np.nan_to_num(dist, nan=0)
        self.dist = dist
        self.axes = axes # list of objects of class Axis
        self.conditional_var = conditional_var # if None this is a joint (or marginal) distribution

        self._validate_distribution()

    def _validate_distribution(self, tol=1e-5):
        # check that probability distributions sum to 1. conditional probability 
        # distributions P(X | Y = y) are allowed to sum to 0 because if P(Y = y) = 0
        # then P(X | Y = y) is nan. 
        if self.conditional_var is not None:
            indexer = ProbabilityDist._index_vars(self.axes)
            self_idx = ProbabilityDist._fetch_idx(self.axes, indexer)
            out_idx = ProbabilityDist._fetch_idx([self.conditional_var], indexer)

            dist_sum = np.einsum(self.dist, self_idx, out_idx)

            if not np.all(np.logical_or(np.isclose(dist_sum, 1, rtol=tol), np.isclose(dist_sum, 0, rtol=tol))):
                raise Exception("help! invalid conditional probability distribution(s) - did not approximately sum to 1\n{}".format(str(dist_sum)))
        else:
            dist_sum = np.sum(self.dist)

            if not np.all(np.isclose(dist_sum, 1, rtol=tol)):
                raise Exception("help! invalid probability distribution(s) - did not approximately sum to 1\n{}".format(str(dist_sum)))
            
        if not np.all(np.logical_or(np.isnan(self.dist), np.logical_and(self.dist <= 1, self.dist >= 0))):
            raise Exception("help! invalid probability distribution(s) - values were outside of [0, 1]")

    @staticmethod
    def _fetch_idx(vars, indexer):
        idx = []
        for x in vars:
            idx.append(indexer[x])
        return(idx)

    @staticmethod
    def _index_vars(vars):
        unique_vars = set(vars)
        var_dict = {var: idx for (idx, var) in enumerate(unique_vars)}
        return(var_dict)

    def condition(self, axis_to_condition):
        """
        Conditions on one of the variables in a joint distribution. 
        e.g. conditioning on y: P(x, y) -> P(x | y) and P(x, y, z) -> P(x, z | y)
        """
        if self.dist.ndim < 2:
            raise Exception("help! cannot generate conditional probabilities from a 1-dimensional distribution")

        if self.conditional_var is not None:
            raise Exception("help! probability distributions cannot condition on more than one varulation")

        # multiply P(x, y) by 1 / P(y) to obtain P(x | y)
        marg_axis = None
        marginal = self
        for axis in self.axes:
            if axis.name != axis_to_condition:
                marginal = marginal.marginalize(axis.name)
            else:
                marg_axis = axis
        if marg_axis is None:
            raise Exception("help! cannot condition on a variable that is not included in this distribution")

        invert_marginal = 1 / marginal.dist

        indexer = ProbabilityDist._index_vars(self.axes)
        self_idx = ProbabilityDist._fetch_idx(self.axes, indexer)
        marg_idx = ProbabilityDist._fetch_idx([marg_axis], indexer)
        out_idx = self_idx

        dist = np.einsum(self.dist, self_idx, invert_marginal, marg_idx, out_idx)

        return(ProbabilityDist(dist, self.axes, conditional_var=marg_axis), marginal)

    def multiply_by_marginal(self, marginal):
        """
        Multiplies a conditional distribution (self) by a marginal distribution
        to yield a joint probability. e.g. P(x | y) * P(y) -> P(x, y) and 
        P(x, z | y) * P(y) -> P(x, y, z)
        """
        if self.conditional_var not in marginal.axes:
            raise Exception("help! marginal distribution does not include conditional variable; cannot generate a joint distribution")

        if self.conditional_var == marginal.conditional_var:
            raise Exception("help! marginal distribution is conditioned on conditional variable; cannot generate a joint distribution")

        indexer = ProbabilityDist._index_vars(self.axes + marginal.axes)
        out_axes = list(indexer.keys())
        self_idx = ProbabilityDist._fetch_idx(self.axes, indexer)
        marg_idx = ProbabilityDist._fetch_idx(marginal.axes, indexer)
        out_idx = ProbabilityDist._fetch_idx(out_axes, indexer)

        dist = np.einsum(self.dist, self_idx, marginal.dist, marg_idx, out_idx)

        return(ProbabilityDist(dist, out_axes, conditional_var=marginal.conditional_var))

    def marginalize(self, axis_to_marginalize):
        """
        Marginalizes out one of the variables in a joint probability. e.g. 
        marginalizing out y: P(x, y) -> P(x) and P(x, y | z) -> P(x | z)
        """
        if self.conditional_var is not None and self.conditional_var.name == axis_to_marginalize:
            raise Exception("help! cannot marginalize out a conditional variable")

        indexer = ProbabilityDist._index_vars(self.axes)
        out_axes = [axis for axis in self.axes if axis.name != axis_to_marginalize]
        self_idx = ProbabilityDist._fetch_idx(self.axes, indexer)
        out_idx = ProbabilityDist._fetch_idx(out_axes, indexer)

        dist = np.einsum(self.dist, self_idx, out_idx)

        return(ProbabilityDist(dist, out_axes, conditional_var=self.conditional_var))

    def joint(self, other):
        """
        Multiply two (marginal) probability distributions to yield a joint probability.
        e.g. P(x) * P(z) -> P(x, z) and P(x | y) * P(z | y) -> P(x, z | y)
        """
        if self.conditional_var != other.conditional_var:
            raise Exception("help! probability distributions not conditioned on the same variables; cannot compute a joint probability")

        indexer = ProbabilityDist._index_vars(self.axes + other.axes)
        out_axes = list(indexer.keys())
        self_idx = ProbabilityDist._fetch_idx(self.axes, indexer)
        other_idx = ProbabilityDist._fetch_idx(other.axes, indexer)
        out_idx = ProbabilityDist._fetch_idx(out_axes, indexer)

        dist = np.einsum(self.dist, self_idx, other.dist, other_idx, out_idx)

        return(ProbabilityDist(dist, out_axes, conditional_var=self.conditional_var))

    def write(self, path):
        axis_str = "\n".join([str(axis) for axis in self.axes])
        if self.conditional_var:
            cond_str = self.conditional_var.name
        else:
            cond_str = "None"
        if self.dist.ndim > 2:
            header = "reshaped {0}d array {1}.\nconditioning on: {2}.\n".format(self.dist.ndim, 
                                                                              self.dist.shape, 
                                                                              cond_str)
            header = header + axis_str
            np.savetxt(path, self.dist.flatten(), header=header)
        else:
            axis_str = "\n".join([str(axis) for axis in self.axes])
            header = "conditioning on: {0}.\n".format(cond_str) + axis_str
            np.savetxt(path, self.dist, header=header)
        return
