import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, combine_pvalues
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gwas')
    parser.add_argument('--control')
    parser.add_argument('--out')
    args = parser.parse_args() 

    gwas = pd.read_csv(args.gwas, sep='\t', usecols=["SNP", "B", "UKB_WBfreq", "CHBfreq", "YRIfreq"])
    gwas = gwas[gwas.UKB_WBfreq >= 0.01] # filters out some lingering variants with frequency < 0.01
    control = pd.read_csv(args.control, sep='\t', usecols=["SNP", "B", "UKB_WBfreq", "CHBfreq", "YRIfreq"])

    UKB_dense_bins = np.array([0, 0.01, 0.015, 0.02, 0.025] + \
        [i / 100 for i in range(3, 20)] + [i / 50 for i in range(10, 50)] + [0.99999, 1]) 
    UKB_sparse_bins = gwas.UKB_WBfreq.quantile([i / 10 for i in range(10)]) 
    B_bins = np.unique(control.B.quantile([i / 19 for i in range(20)]))
    YRI_bins = np.array([i / 216 for i in range(217)])
    CHB_bins = np.array([i / 206 for i in range(207)])
    
    # bin frequencies and B-values for gwas and control variants
    gwas["UKB_dense_bin"] = np.digitize(gwas.UKB_WBfreq, bins=UKB_dense_bins)
    gwas["UKB_sparse_bin"] = np.digitize(gwas.UKB_WBfreq, bins=UKB_sparse_bins)
    gwas["B_bin"] = np.digitize(gwas.B, bins=B_bins)
    gwas["YRI_bin"] = np.digitize(gwas.YRIfreq, bins=YRI_bins) - 1
    gwas["CHB_bin"] = np.digitize(gwas.CHBfreq, bins=CHB_bins) - 1
    control["UKB_dense_bin"] = np.digitize(control.UKB_WBfreq, bins=UKB_dense_bins)
    control["UKB_sparse_bin"] = np.digitize(control.UKB_WBfreq, bins=UKB_sparse_bins)
    control["B_bin"] = np.digitize(control.B, bins=B_bins)
    control["YRI_bin"] = np.digitize(control.YRIfreq, bins=YRI_bins) - 1
    control["CHB_bin"] = np.digitize(control.CHBfreq, bins=CHB_bins) - 1
    
    # count number of control variants in each combination of (UKB WB frequency and B value) bins    
    matched_counts = control.groupby(["UKB_dense_bin", "B_bin"]).count().rename({"UKB_sparse_bin": "n_matched"}, axis=1)["n_matched"]
    
    # drop gwas variants with insufficient matched variants
    enough_matched_counts = matched_counts[matched_counts >= 500]
    gwas = pd.merge(enough_matched_counts, gwas, left_index=True, right_on=["UKB_dense_bin", "B_bin"])
    gwas.n_matched.to_csv(args.out + "_matched_counts.txt", sep='\t', index=False)
    
    print("starting merge...")
    # match variants via an incredibly memory-intensive merge
    # (there are probably better ways to do this but alas)
    matched = pd.merge(control, gwas[["UKB_dense_bin", "B_bin", "SNP", "YRIfreq"]].rename({"SNP": "gwas_SNP", "YRIfreq": "gwas_YRI"}, axis=1), 
                       left_on=["UKB_dense_bin", "B_bin"], right_on=["UKB_dense_bin", "B_bin"])
    matched = matched[["CHBfreq", "YRIfreq", "UKB_sparse_bin", "UKB_WBfreq", "YRI_bin", "CHB_bin"]]
    print("finished merge...")
    
    # generate CFS for gwas and matched variants
    gwas_YRI = gwas.groupby(["UKB_sparse_bin", "YRI_bin"]).count().reset_index()[["UKB_sparse_bin", "YRI_bin", "UKB_WBfreq"]].rename({"UKB_WBfreq": "count"})
    gwas_CHB = gwas.groupby(["UKB_sparse_bin", "CHB_bin"]).count().reset_index()[["UKB_sparse_bin", "CHB_bin", "UKB_WBfreq"]].rename({"UKB_WBfreq": "count"})
    matched_YRI = matched.groupby(["UKB_sparse_bin", "YRI_bin"]).count().reset_index()[["UKB_sparse_bin", "YRI_bin", "UKB_WBfreq"]].rename({"UKB_WBfreq": "count"})
    matched_CHB = matched.groupby(["UKB_sparse_bin", "CHB_bin"]).count().reset_index()[["UKB_sparse_bin", "CHB_bin", "UKB_WBfreq"]].rename({"UKB_WBfreq": "count"})
    gwas_YRI.to_csv(args.out + "_gwas_YRI_cfs.txt", sep='\t', index=False)
    gwas_CHB.to_csv(args.out + "_gwas_CHB_cfs.txt", sep='\t', index=False)
    matched_YRI.to_csv(args.out + "_matched_YRI_cfs.txt", sep='\t', index=False)
    matched_CHB.to_csv(args.out + "_matched_CHB_cfs.txt", sep='\t', index=False)

    # compute heterozygosity
    gwas["CHBhet"] =  2 * gwas.CHBfreq * (1 - gwas.CHBfreq)
    gwas["YRIhet"] = 2 * gwas.YRIfreq * (1 - gwas.YRIfreq)
    matched["CHBhet"] = 2 * matched.CHBfreq * (1 - matched.CHBfreq)
    matched["YRIhet"] = 2 * matched.YRIfreq * (1 - matched.YRIfreq)
    
    # compute pvalues in each bin and aggregate
    pvals = pd.DataFrame({"bin": gwas.UKB_sparse_bin.unique()})
    for pop in ["CHB", "YRI"]:
        two_sided_tests = []
        one_sided_tests = []
        heterozygosity_pvals = []
        for i in gwas.UKB_sparse_bin.unique():
            x = ttest_ind(gwas[gwas.UKB_sparse_bin == i][pop + "freq"], matched[matched.UKB_sparse_bin == i][pop + "freq"], equal_var=False)
            two_sided_tests.append(x)
            y = ttest_ind(gwas[gwas.UKB_sparse_bin == i][pop + "freq"], matched[matched.UKB_sparse_bin == i][pop + "freq"], equal_var=False, alternative="less")
            one_sided_tests.append(x)
            z = ttest_ind(gwas[gwas.UKB_sparse_bin == i][pop + "het"], matched[matched.UKB_sparse_bin == i][pop + "het"], equal_var=False)
            heterozygosity_pvals.append(z.pvalue)
        two_sided_pvals, two_sided_stats = [x.pvalue for x in two_sided_tests], [x.statistic for x in two_sided_tests]
        two_sided_combined_pval = combine_pvalues(two_sided_pvals, method="fisher")
        one_sided_pvals, one_sided_stats = [x.pvalue for x in one_sided_tests], [x.statistic for x in one_sided_tests]
        one_sided_combined_pval = combine_pvalues(one_sided_pvals, method="fisher")
        conservative_one_sided_combined_pval = combine_pvalues(one_sided_pvals[:2], method="fisher")
        heterozygosity_combined_pval = combine_pvalues(heterozygosity_pvals, method="fisher")
        pvals[pop + "_two_sided_statistic"] = two_sided_stats
        pvals[pop + "_two_sided_pvalue"] = two_sided_pvals
        pvals[pop + "_two_sided_meta_pvalue"] = two_sided_combined_pval[1]
        pvals[pop + "_one_sided_statistic"] = one_sided_stats
        pvals[pop + "_one_sided_pvalue"] = one_sided_pvals
        pvals[pop + "_one_sided_meta_pvalue"] = one_sided_combined_pval[1]
        pvals[pop + "_conservative_one_sided_meta_pvalue"] = conservative_one_sided_combined_pval[1]
        pvals[pop + "_heterozygosity_meta_pvalue"] = heterozygosity_combined_pval[1]
    pvals.to_csv(args.out + "_gwas_pvalues.txt", sep='\t', index=False)

    # compute mean and bootstrap over bins to compute standard errors
    matched_summary = pd.concat([matched.groupby("UKB_sparse_bin")[["UKB_WBfreq"]].median().add_suffix("_median"), 
                                 matched.groupby("UKB_sparse_bin")[["CHBfreq", "YRIfreq", "CHBhet", "YRIhet"]].mean().add_suffix("_mean")], axis=1)
    matched_samples = []
    for i in range(100):
        print(i)
        sample = matched.groupby("UKB_sparse_bin")[["YRIfreq", "CHBfreq", "YRIhet", "CHBhet"]].agg(lambda x: x.sample(frac = 1.0, replace = True).mean())
        matched_samples.append(sample) 
    matched_samples = pd.concat(matched_samples).groupby("UKB_sparse_bin").agg([('lower', lambda x: x.quantile(0.025)), ('upper', lambda x: x.quantile(0.975))])
    matched_samples.columns = matched_samples.columns.map('_'.join)
    matched_summary = pd.concat([matched_summary, matched_samples], axis=1)
    matched_summary.to_csv(args.out + "_matched_summary.txt", sep='\t', index=False)

    gwas_summary = pd.concat([gwas.groupby("UKB_sparse_bin")[["UKB_WBfreq"]].median().add_suffix("_median"), 
                              gwas.groupby("UKB_sparse_bin")[["CHBfreq", "YRIfreq", "CHBhet", "YRIhet"]].mean().add_suffix("_mean")], axis=1)
    gwas_samples = []
    for _ in range(100):
        sample = gwas.groupby("UKB_sparse_bin")[["YRIfreq", "CHBfreq", "YRIhet", "CHBhet"]].agg(lambda x: x.sample(frac = 1.0, replace = True).mean())
        gwas_samples.append(sample) 
    gwas_samples = pd.concat(gwas_samples).groupby("UKB_sparse_bin").agg([('lower', lambda x: x.quantile(0.025)), ('upper', lambda x: x.quantile(0.975))])
    gwas_samples.columns = gwas_samples.columns.map('_'.join)
    gwas_summary = pd.concat([gwas_summary, gwas_samples], axis=1)
    gwas_summary.to_csv(args.out + "_gwas_summary.txt", sep='\t', index=False)

if __name__ == "__main__":
    main()
