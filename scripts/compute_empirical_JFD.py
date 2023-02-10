import numpy as np
from jfd import EmpiricalDist 
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gwas', action="store_true")
    parser.add_argument('--control', action="store_true")
    parser.add_argument('--data')
    parser.add_argument('--condition', nargs='+')
    parser.add_argument('--control_jfd', nargs='*')
    parser.add_argument('--B_bounds', default=None)
    parser.add_argument('--pop')
    parser.add_argument('--out_counts')
    parser.add_argument('--out_jfd')
    parser.add_argument('--out_null')
    args = parser.parse_args()
    
    print(args.condition)

    if args.B_bounds is not None:
        B_bin_bounds = np.loadtxt(args.B_bounds)
    bounds = {}
    for var in args.condition:
        if var == "B":
            bounds[var] = B_bin_bounds
        elif var == EmpiricalDist.study_pop:
            bounds[var] = EmpiricalDist.study_pop_bounds
        else:
            bounds[var] = EmpiricalDist.minor_pop_bounds

    empirical = EmpiricalDist(args.condition, bounds)
    empirical.add_observations(args.data)
    empirical.initialize()

    if args.gwas and not args.control:
        for fp in args.control_jfd:
            control_dist = np.load(fp) 
            null = empirical.compute_null(control_dist)

            tbl = np.concatenate([empirical.data, null.T], axis=1)

            pop = fp[-7:-4]
            out_fp = args.out_null + pop + ".txt"
            with open(out_fp, 'a') as f:
                f.write("\t".join(empirical.header + [str(i / 100) for i in range(101)]))
                f.write("\n")
                np.savetxt(f, tbl, delimiter='\t')

        _, dist = empirical.compute_study_conditional_distribution(args.pop + "freq")
        np.save(args.out_jfd + args.pop + ".npy", dist)

#         for pop in empirical.pops:
#             if pop != empirical.study_pop:
#                 _, dist = empirical.compute_study_conditional_distribution(pop)
#                 np.save(args.out_jfd + pop[:-4] + ".npy", dist)                

    elif args.control and not args.gwas:
        counts, dist = empirical.compute_study_conditional_distribution(args.pop + "freq")
        counts = np.sum(counts, axis=0)
        np.save(args.out_counts, counts)
        np.save(args.out_jfd + args.pop + ".npy", dist)
#         stored_counts = False
#         for pop in empirical.pops:
#             if pop != empirical.study_pop:
#                 counts, dist = empirical.compute_study_conditional_distribution(pop)
#                 if not stored_counts: 
#                     counts = np.sum(counts, axis=0)
#                     np.save(args.out_counts, counts)
#                     stored_counts = True
#                 np.save(args.out_jfd + pop[:-4] + ".npy", dist)                
    
if __name__ == "__main__":
    main()