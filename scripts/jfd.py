import numpy as np

class JointFrequencyDist:
    bin_bounds = np.linspace(0, 1.01, 102)

    def __init__(self) -> None:
        self.data = np.empty((0,0))
        self.header = []
        # self.sample_bins = None
        # self.bin_counts = None
        self.initialized = False

    def add_observations(self, filepath=None, arr=None, header=None):
        if filepath is not None:
            with open(filepath, 'r') as f:
                header = f.readline().strip().replace('"', '').split('\t')
            new_data = np.loadtxt(filepath, skiprows = 1)
        elif arr is not None and header is not None:
            new_data = arr

        # if "B" in header: # Need to reorder data columns so that B-value comes last
        #     B_idx = header.index("B")
        #     order = list(range(len(header)))
        #     reorder = order[:B_idx] + order[B_idx + 1:] + [order[B_idx]]
        #     header = [header[i] for i in reorder]
        #     new_data = new_data[:, reorder]

        if not self.header:
            self.header = header
        elif self.header != header:
            raise Exception("population mismatch between existing {0} and new {1} data".format(self.header, header))

        if not np.any(self.data):
            self.data = new_data
        else:
            self.data = np.concatenate([self.data, new_data], axis=0)

    def initialize(self):
        self.n = self.data.shape[0]
        self.initialized = True

    def idx(self, val):
        return(self.header.index(val))

    def summary(self):
        if not self.initialized:
            self.initialize()
        print("Num. observations: ", self.n)
        if self.header is None:
            print("Populations: None")
        else:
            print("Populations: ", ", ".join(self.header))

    @staticmethod
    def normalize(arr):
        if len(arr.shape) == 1:
            probs = arr / np.sum(arr)
        elif len(arr.shape) > 1 and len(arr.shape) < 4:
            probs = arr / np.sum(arr, axis=0)
        else:
            raise Exception("matrix dimensions incorrect for JFD.normalize()")
        return(probs)

    ### DEPRECATED (not used afaik) ###
    def bin_frequencies(self, pop):
        """Assign each observation to a frequency bin based on frequency in input 'pop'.
        Return array of bins and, if 'pop' is focal pop, set self.sample_bins equal to this array."""
        idx = self.pops.index(pop)
        bins = np.digitize(self.data[:,idx], bins = self.freq_bin_bounds) - 1
        counts, _ = np.histogram(self.data[:,idx], bins = self.freq_bin_bounds)
        if np.any(bins >= self.n_bins) or np.any(bins < 0):
            raise Exception("frequencies out of bin boundary range")
        if pop == self.focal_pop:
            self.sample_bins = bins
            self.bin_counts = counts
        return(bins, counts)

    ### DEPRECATED (not used afaik) ###
    def fetch_bin_probabilities(self):
        idx = self.pops.index(self.focal_pop)
        bin_probs, _ = np.histogram(self.data[:,idx], bins = self.freq_bin_bounds, density = True)
        return(bin_probs)

    ### DEPRECATED (not used afaik) ###
    @staticmethod
    def reindex_vals(bin_vals, sample_bins, starting_idx = 0):
        # note to self - starting_idx must explicitly be passed in because we can't just 
        # assume the min value of sample_bins suffices as the starting_idx, that 
        # assumption will often not be true
        sample_vals = bin_vals[sample_bins - starting_idx]
        return(sample_vals)

class EmpiricalDist(JointFrequencyDist):
    study_pop = "UKB_WBfreq"
    study_pop_bounds = JointFrequencyDist.bin_bounds
    minor_pop_bounds = np.linspace(0, 1.05, 22)

    def __init__(self, cond_var, bounds_dict):
        JointFrequencyDist.__init__(self)
        self.conditioning_var = cond_var
        self.bounds = bounds_dict
        # self.study_conditional_dist = None ####### P(x_k | x_S, B) [k x m x m x b matrix; 0 if k == S]

    def initialize(self):
        JointFrequencyDist.initialize(self)
        # self.study_conditional_dist = np.full((self.k, self.m, self.m), 0.0)
        # self.B_study_conditional_dist = np.full((self.k, self.m, self.m, self.b), 0.0)

    def summary(self):
        JointFrequencyDist.summary(self)
        print("Conditioning variables: ", ", ".join(self.conditioning_var))

    def compute_study_conditional_distribution(self, pop):
        """Compute empirical probability distribution of `pop` conditional on study population, i.e.
           P(x_k | x_S) where k != S. Can optionally condition on additional variables passed in *args;
           *args should contain a list of tuples where the first element is the name of the variable 
           to condition on (in self.header) and the second element is the bin bounds to use."""

        if not self.initialized:
            self.initialize()

        data_list = [self.data[:, self.idx(pop)]]
        bounds_list = [self.study_pop_bounds]
        for var in self.conditioning_var:
            data_list.append(self.data[:, self.idx(var)])
            bounds_list.append(self.bounds[var])

        counts, _ = np.histogramdd(np.array(data_list).T, bins=bounds_list)
        dist = self.normalize(counts) # P(x_k | x_S, B)
        # self.B_study_conditional_dist[self.idx(pop)] = dist
        # self.study_conditional_dist[self.idx(pop)] = np.sum(dist, axis=2)
        return(counts, dist)

    def bin_samples(self, var):
        """Assign each observation to a bin based on varue of 'var'. Return array of bins."""
        idx = self.header.index(var)
        bounds = self.bounds[var]
        bins = np.digitize(self.data[:,idx], bins = bounds) - 1
        if np.any(bins >= len(bounds) - 1) or np.any(bins < 0):
            raise Exception("varues out of bin boundary range")
        return(bins)

    def compute_null(self, control_dist):
        if len(self.conditioning_var) == 1:
            bins = self.bin_samples(self.conditioning_var[0])
            null = control_dist[:, bins]
        elif len(self.conditioning_var) == 2:
            bin_1 = self.bin_samples(self.conditioning_var[0])
            bin_2 = self.bin_samples(self.conditioning_var[1])
            null = control_dist[:, bin_1, bin_2]
        elif len(self.conditioning_var) == 3:
            bin_1 = self.bin_samples(self.conditioning_var[0])
            bin_2 = self.bin_samples(self.conditioning_var[1])
            bin_3 = self.bin_samples(self.conditioning_var[2])
            null = control_dist[:, bin_1, bin_2, bin_3]
        self.null = null
        return(null)

    ### DEPRECATED (not used afaik) ###
    def set_B_bin_bounds(self, fp):
        self.B_bin_bounds = np.loadtxt(fp)
        self.b = len(self.B_bin_bounds) - 1

    ### DEPRECATED (not used afaik) ###
    def compute_likelihood(self, prob_dist, pop):
        if not np.any(self.study_conditional_counts[self.idx(pop)]):
            self.compute_study_conditional_distribution(pop)

        prob_dist[np.where(prob_dist == 0)] = np.finfo(float).eps # converts 0 to epsilon to ensure log transformation results in real-valued numbers
        log_prob = np.log(prob_dist)
        composite_log_likelihood = np.sum(np.multiply(self.study_conditional_counts[self.idx(pop)], log_prob))
        return(composite_log_likelihood)

    ### DEPRECATED (not used afaik) ###
    def fetch_quantiles(self, prob_dist, pop):
        cdf = np.cumsum(prob_dist, axis=0)[:,10:]
        counts = self.study_conditional_counts[self.idx(pop)][:,10:]
        quantiles = np.repeat(cdf.flatten(order='F'), counts.flatten(order='F'))
        return(quantiles)

class SimulatedDist(JointFrequencyDist):
    freq_bin_bounds = JointFrequencyDist.bin_bounds
    study_pop = "CEUfreq"
    ancestral_pop = "ancestral_freq"
    differential = 0.01

    def __init__(self):
        JointFrequencyDist.__init__(self)
        self.ancestral_dist = None ############### P(x_A) [m x 1 array]
        self.ancestral_conditional_dist = None ### P(x_k | x_A) [k x m x m matrix; 0 if k == A]
        self.marginal_study_dist = None ########## P(x_S) [m x 1 array]
        self.study_conditional_dist = None ####### P(x_k | x_S) [k x m x m matrix; 0 if k == S]

    def initialize(self):
        JointFrequencyDist.initialize(self)
        self.k = len(self.header)
        self.m = len(self.freq_bin_bounds) - 1
        self.ancestral_conditional_dist = np.full((self.k, self.m, self.m), 0.0)
        self.study_conditional_dist = np.full((self.k, self.m, self.m), 0.0)
        self.ancestral_dist = np.full((self.m), 0.0)
        self.marginal_study_dist = np.full((self.m), 0.0)

    def summary(self):
        JointFrequencyDist.summary(self)
        print("Study population: ", self.study_pop)

    def set_ancestral_distribution(self, arr):
        self.ancestral_dist = arr
        return(arr)

    def resample_data(self, pop, n):
        sample_counts = np.random.binomial(n, self.data[:,self.idx(pop)])
        sample_freqs = sample_counts / n
        self.data[:,self.idx(pop)] = sample_freqs
        return(self.data[:,self.idx(pop)])

    def compute_ancestral_conditional_distribution(self, pop):
        """Compute probability distribution of `pop` conditional on ancestral population, i.e.
           P(x_k | x_A) where k != A."""

        if not self.initialized:
            self.initialize()

        y = self.data[:, self.idx(self.ancestral_pop)]
        x = self.data[:, self.idx(pop)]
        counts, x_, y_ = np.histogram2d(x, y, bins=self.freq_bin_bounds)
        dist = self.normalize(counts) # P(x_k | x_A)
        self.ancestral_conditional_dist[self.idx(pop)] = dist
        return(dist)

    def compute_marginal_study_distribution(self):
        """Compute marginal distribution of study population, i.e. P(x_S)."""

        if not self.initialized:
            self.initialize()
        if not np.any(self.ancestral_conditional_dist[self.idx(self.study_pop)]):
            self.compute_ancestral_conditional_distribution(self.study_pop)
        if not np.any(self.ancestral_dist):
            raise Exception("no ancestral distribution provided")

        study_dist = self.ancestral_conditional_dist[self.idx(self.study_pop)] # P(x_S | x_A)
        study_dist = np.nan_to_num(study_dist) # replace nan with 0 for matrix math
        dist = np.matmul(study_dist, self.ancestral_dist) # P(x_S)
        self.marginal_study_dist = dist
        return(dist)

    def compute_study_conditional_distribution(self, pop):
        """Compute probability distribution of `pop` conditional on study population, i.e.
           P(x_k | x_S) where k != S. Note that the math is different for k == A vs k != A."""

        if not self.initialized:
            self.initialize()

        if pop == self.ancestral_pop:
            if not np.any(self.marginal_study_dist):
                self.compute_marginal_study_distribution()
            if not np.any(self.ancestral_dist):
                raise Exception("no ancestral distribution provided")

            study_dist = self.ancestral_conditional_dist[self.idx(self.study_pop)] # P(x_S | x_A)
            study_dist = np.nan_to_num(study_dist) # replace nan with 0 for matrix math
            dist = np.divide(np.multiply(study_dist, self.ancestral_dist).T,
                             self.marginal_study_dist) # P(x_A | x_S)
        else:
            if not np.any(self.study_conditional_dist[self.idx(self.ancestral_pop)]):
                self.compute_study_conditional_distribution(self.ancestral_pop)
            if not np.any(self.ancestral_conditional_dist[self.idx(pop)]):
                self.compute_ancestral_conditional_distribution(pop)

            pop_dist = self.ancestral_conditional_dist[self.idx(pop)] # P(x_k | x_A)
            pop_dist = np.nan_to_num(pop_dist) # replace nan with 0 for matrix math
            ancestral_dist = self.study_conditional_dist[self.idx(self.ancestral_pop)] # P(x_A | x_S)
            ancestral_dist = np.nan_to_num(ancestral_dist) # replace nan with 0 for matrix math
            dist = np.matmul(pop_dist, ancestral_dist) # P(x_k | x_S)
            
        self.study_conditional_dist[self.idx(pop)] = dist
        return(dist)
