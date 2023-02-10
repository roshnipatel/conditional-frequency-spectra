import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dist')
    parser.add_argument('--out')
    args = parser.parse_args()

    null = pd.read_csv(args.dist, sep='\t')

    probs = null.iloc[:,-101:]
    freqs = [i / 100 for i in range(101)]
    null["mean"] = (probs * freqs).sum(axis=1)
    null["variance"] = probs.apply(lambda row: np.sum(row * (freqs - null.iloc[row.name]["mean"]) ** 2), axis=1)

    null[["sample_" + str(i) for i in range(10)]] = \
        probs.apply(lambda row: np.random.choice(freqs, 10, p=row), axis=1, result_type='expand')
    pop = args.dist[-7:-4]
    null["gwas_z"] = (null[pop + "freq"] - null["mean"]) / np.sqrt(null["variance"])
    null["gwas_p"] = norm.sf(abs(null["gwas_z"])) * 2

    null[["p_" + str(i) for i in range(10)]] = \
        null.apply(lambda row: norm.sf(abs((row[["sample_" + str(i) for i in range(10)]] - row["mean"]) / np.sqrt(row["variance"]))) * 2, 
                axis=1, result_type='expand')

    null.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()