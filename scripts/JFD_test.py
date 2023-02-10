import unittest
import numpy as np
from frequency_distribution import JointFrequencyDist, SimulatedDist, EmpiricalDist

class TestSimJFD(unittest.TestCase):
    @staticmethod
    def jfd_example():
        jfd = SimulatedDist()
        jfd.add_observations(arr = np.array([[0.0101, 0, 0],
                                             [0.501, 0.301, 0.601],
                                             [0.801, 1, 0.901],
                                             [0.701, 0.301, 0.201]]),
                             pops = ["CEUfreq", "ancestral_freq", "CHBfreq"])
        return(jfd)

    @staticmethod
    def anc_example():
        anc_dist = np.full((101), 0.0)
        anc_dist[[0, 30, 100]] = 1/3
        return(anc_dist)

    @staticmethod
    def empirical_example():
        jfd = EmpiricalDist()
        jfd.add_observations(arr = np.array([[0.0101, 0],
                                             [0.0101, 0],
                                             [0.501, 0.201],
                                             [0.501, 0.601],
                                             [0.701, 0.201],
                                             [0.701, 0.201],
                                             [0.801, 0.901],
                                             [0.801, 0.901],
                                             [0.801, 0.901],
                                             [0.801, 0.901]]),
                             pops = ["CEUfreq", "CHBfreq"])
        return(jfd)

    @staticmethod
    def jfd_example2():
        jfd = SimulatedDist()
        jfd.add_observations(arr = np.array([[0.101, 0.501, 0.501],
                                             [0.401, 0.501, 0.401],
                                             [0.301, 0.301, 0.401],
                                             [0.401, 0.301, 0.201]]),
                             pops = ["ancestral_freq", "CEUfreq", "CHBfreq"])
        return(jfd)

    @staticmethod
    def anc_example2():
        anc_dist = np.full((101), 0.0)
        anc_dist[10] = 0.7
        anc_dist[40] = 0.2
        anc_dist[30] = 0.1
        return(anc_dist)

    def compare_arrays(self, arr1, arr2):
        def check_equality(flat1, flat2):
            for i in range(len(flat1)):
                self.assertAlmostEqual(flat1[i], flat2[i])

        self.assertEqual(arr1.shape, arr2.shape)
        if not (np.all(arr1) and np.all(arr2)):
            nan1 = np.isnan(arr1)
            nan2 = np.isnan(arr2)
            check_equality(nan1.flatten(), nan2.flatten())
            check_equality(arr1[~nan1].flatten(), arr2[~nan2].flatten())
        else:
            check_equality(arr1, arr2)

class TestEmpirical(TestSimJFD):
    def test_pop_study(self):
        jfd = self.empirical_example()
        jfd.initialize()
        expected = np.full((jfd.m, jfd.m), 0.0)
        expected[0, 1] = 2
        expected[[60, 20], 50] = 1
        expected[20, 70] = 2
        expected[90, 80] = 4

        jfd.compute_study_conditional_distribution("CHBfreq")
        output = jfd.study_conditional_counts[jfd.idx("CHBfreq")]

        self.compare_arrays(output, expected)

    def test_quantiles(self):
        jfd = self.jfd_example()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example())
        cdf = jfd.compute_study_conditional_distribution("CHBfreq")

        empirical = self.empirical_example()
        empirical.compute_study_conditional_distribution("CHBfreq")

        expected = np.array([1, 1, 0.5, 1, 0.5, 0.5, 1, 1, 1, 1])
        output = empirical.fetch_quantiles(cdf, "CHBfreq")

        self.compare_arrays(output, expected)

class TestBasicFns(TestSimJFD):
    def test_initialization(self):
        jfd = self.jfd_example()
        jfd.initialize()
        self.assertEqual(jfd.n, 4)
        self.assertEqual(jfd.k, 3)
        self.assertEqual(jfd.m, 101)

    def test_resample(self):
        jfd = self.jfd_example()
        CEU_initial = np.copy(jfd.data[:,jfd.idx("CEUfreq")])
        jfd.initialize()
        CEU_sampled = jfd.resample_data("CEUfreq", 10**16)
        self.compare_arrays(CEU_initial, CEU_sampled)

class TestExample1(TestSimJFD):
    def test_anc_cond(self):
        jfd = self.jfd_example()
        jfd.initialize()
        output = jfd.compute_ancestral_conditional_distribution("CEUfreq")
        expected = np.empty((jfd.m, jfd.m))
        expected[:] = np.nan
        expected[:, [0, 30, 100]] = 0
        expected[1, 0] = 1
        expected[[50, 70], 30] = 0.5
        expected[80, 100] = 1
        self.compare_arrays(output, expected)

    def test_marg_study(self):
        jfd = self.jfd_example()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example())
        output = jfd.compute_marginal_study_distribution()
        expected = np.full((jfd.m), 0.0)
        expected[[1, 80]] = 1/3
        expected[[70, 50]] = 1/6
        self.compare_arrays(output, expected)

    def test_anc_study(self):
        jfd = self.jfd_example()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example())
        output = jfd.compute_study_conditional_distribution(jfd.ancestral_pop)
        expected = np.empty((jfd.m, jfd.m))
        expected[:] = np.nan
        expected[:, [1, 50, 70, 80]] = 0
        expected[0, 1] = 1
        expected[30, [50, 70]] = 1
        expected[100, 80] = 1
        self.compare_arrays(output, expected)

    def test_pop_study(self):
        jfd = self.jfd_example()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example())
        output = jfd.compute_study_conditional_distribution("CHBfreq")
        expected = np.full((jfd.m, jfd.m), 0.0)
        expected[0, 1] = 1
        expected[[60, 20], 50] = 0.5
        expected[[60, 20], 70] = 0.5
        expected[90, 80] = 1
        self.compare_arrays(output, expected)

class TestExample2(TestSimJFD):
    def test_anc_cond(self):
        jfd = self.jfd_example2()
        jfd.initialize()
        output = jfd.compute_ancestral_conditional_distribution("CEUfreq")
        expected = np.empty((jfd.m, jfd.m))
        expected[:] = np.nan
        expected[:, [10, 30, 40]] = 0
        expected[50, 10] = 1
        expected[[30, 50], 40] = 0.5
        expected[30, 30] = 1
        self.compare_arrays(output, expected)
    
    def test_marg_study(self):
        jfd = self.jfd_example2()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example2())
        output = jfd.compute_marginal_study_distribution()
        expected = np.full((jfd.m), 0.0)
        expected[50] = 0.8
        expected[30] = 0.2
        self.compare_arrays(output, expected)

    def test_anc_study(self):
        jfd = self.jfd_example2()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example2())
        output = jfd.compute_study_conditional_distribution(jfd.ancestral_pop)
        expected = np.empty((jfd.m, jfd.m))
        expected[:] = np.nan
        expected[:, [30, 50]] = 0
        expected[10, 50] = 7/8
        expected[40, 50] = 1/8
        expected[[30, 40], 30] = 0.5
        self.compare_arrays(output, expected)

    def test_pop_study(self):
        jfd = self.jfd_example2()
        jfd.initialize()
        jfd.set_ancestral_distribution(self.anc_example2())
        output = jfd.compute_study_conditional_distribution("CHBfreq")
        expected = np.full((jfd.m, jfd.m), 0.0)
        expected[50, 50] = 7/8
        expected[[40, 20], 50] = 1/16
        expected[40, 30] = 3/4
        expected[20, 30] = 1/4
        self.compare_arrays(output, expected)

if __name__ == '__main__':
    unittest.main()