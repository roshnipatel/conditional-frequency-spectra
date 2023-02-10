import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gwas')
    parser.add_argument('--snps')
    parser.add_argument('--out_gwas')
    parser.add_argument('--out_control')
    args = parser.parse_args()
    
    gwas = pd.read_csv(args.gwas, sep=',').drop_duplicates()
    snps = pd.read_csv(args.snps, sep='\t', names=["rsID", "ancestral", "SNP", "alt_freq", "B"])
    
    gwas = pd.merge(gwas, snps, on="SNP")
    gwas.to_csv(args.out_gwas, sep='\t', columns=["SNP", "rsID", "alt_freq", "ancestral", "B"], index=False)
    
    snps = snps[~snps["SNP"].isin(gwas["SNP"])]
    snps.to_csv(args.out_control, sep='\t', columns=["SNP", "rsID", "alt_freq", "ancestral", "B"], index=False)

if __name__ == "__main__":
    main()