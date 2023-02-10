import pandas as pd
import argparse

def merge_freq_file(anc, freq, pop):
    vals_to_remove = freq.apply(lambda row: '<' in row.a1 or '<' in row.a2, axis=1) # cleaning up data
    freq = freq[~vals_to_remove]
    freq[["allele1", "freq1"]] = freq.apply(lambda row: row.a1.split(':'), axis=1, result_type='expand')
    freq[["allele2", "freq2"]] = freq.apply(lambda row: row.a2.split(':'), axis=1, result_type='expand')
    merged = pd.merge(anc, freq, on=["chrom", "pos"])
    a1_ancestral = merged.loc[merged.allele1 == merged.ancestral,:]
    a1_ancestral[pop + "freq"] = a1_ancestral["freq2"]
    a2_ancestral = merged.loc[merged.allele2 == merged.ancestral,:]
    a2_ancestral[pop + "freq"] = a2_ancestral["freq1"]
    merged = pd.concat([a1_ancestral, a2_ancestral]).drop(["n_alleles", "n_chr", "a1", "a2", "allele1", "allele2", "freq1", "freq2"], axis=1)
    return(merged)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ancestral_file')
    parser.add_argument('--frequency_files', nargs='+')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    anc_table = pd.read_csv(args.ancestral_file, sep='\t', names=["SNP", "rsID", "alt_freq", "ancestral", "B"])
    anc_table[["chrom", "pos", "ref", "alt"]] = anc_table.apply(lambda row: row.SNP.split(':'), axis=1, result_type="expand")
    anc_table["chrom"] = pd.to_numeric(anc_table["chrom"])
    anc_table["pos"] = pd.to_numeric(anc_table["pos"])
    ref_ancestral = anc_table.loc[anc_table.ref == anc_table.ancestral,:]
    ref_ancestral["UKB_WBfreq"] = ref_ancestral["alt_freq"]
    alt_ancestral = anc_table.loc[anc_table.alt == anc_table.ancestral,:]
    alt_ancestral["UKB_WBfreq"] = 1 - alt_ancestral["alt_freq"]
    anc_table = pd.concat([ref_ancestral, alt_ancestral])
    
    for f in args.frequency_files:
        freq_table = pd.read_csv(f, sep='\t', header=0, names=["chrom", "pos", "n_alleles", "n_chr", "a1", "a2"])
        pop_idx = f.find("pop") + 3
        pop = f[pop_idx:pop_idx + 3]
        anc_table = merge_freq_file(anc_table, freq_table, pop)
    
    anc_table.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()
