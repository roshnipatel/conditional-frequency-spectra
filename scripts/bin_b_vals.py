import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data')
    parser.add_argument('--out')
    args = parser.parse_args()

    with open(args.data, 'r') as f:
        header = f.readline().strip().replace('"', '').split('\t')
    
    B_idx = header.index("B")
    B_vals = np.loadtxt(args.data, skiprows=1, usecols=B_idx)
    quants = np.unique(np.quantile(B_vals, np.linspace(0, 1, 20))) # 15 bins
    np.savetxt(args.out, quants)

if __name__ == "__main__":
    main()