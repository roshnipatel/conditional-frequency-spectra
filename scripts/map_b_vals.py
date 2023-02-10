import numpy as np
import csv
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bvals')
    parser.add_argument('--snps')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    b_val = np.loadtxt(args.bvals, usecols=0)
    b_range = np.loadtxt(args.bvals, usecols=1)
    b_pos = np.cumsum(b_range)
    
    with open(args.snps, 'r') as f:
        features = f.readline().strip().replace('"', '').split('\t')
    snps = np.loadtxt(args.snps, usecols=features.index("position"), skiprows=1) # must be sorted by position!
    
    snp_b_vals = []
    b_pos_ptr = 0
    for snp_ptr in range(snps.shape[0]):
        while snps[snp_ptr] > b_pos[b_pos_ptr]:
            b_pos_ptr += 1
        snp_b_vals.append(b_val[b_pos_ptr])
    print("snp_b_vals: ", snp_b_vals)
    
    read_file, write_file = open(args.snps, 'r'), open(args.out, 'w')
    reader = csv.reader(read_file, delimiter='\t')
    writer = csv.writer(write_file, delimiter='\t')
    headers = reader.__next__()
    headers.append("B")
    writer.writerow(headers)
    ptr = 0
    for row in reader:
        row.append(snp_b_vals[ptr])
        writer.writerow(row)
        ptr += 1
    read_file.close()
    write_file.close()

if __name__ == "__main__":
    main()