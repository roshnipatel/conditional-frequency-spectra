import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snps')
    parser.add_argument('--table')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    snps = pd.read_csv(args.snps, sep='\t')
    rsIDs = set(snps["rsID"])

    with open(args.table, 'r') as fr:
        with open(args.out, 'w') as fw:
            for line in fr:
                curr_rsID = line.strip().split('\t')[1]
                if curr_rsID in rsIDs:
                    fw.write(line)

if __name__ == "__main__":
    main()