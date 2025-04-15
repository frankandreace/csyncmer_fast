#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from sys import argv


if __name__ == "__main__" :
    
    file = argv[1]
    outfile = argv[2]

    df = pd.read_csv(file,delimiter="\t")

    all_columns = ["HASHING","NAIVE","DEQUE","SYNG_ORIGINAL","RESCAN","RESCAN_CA_BRANCHLESS","RESCAN_CA","RESCAN_CA_ITERATOR"]

    missing_columns = [col for col in all_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"The following columns are not in the DataFrame: {missing_columns}")

    syncmer_columns = ["NAIVE","DEQUE","SYNG_ORIGINAL","RESCAN","RESCAN_CA_BRANCHLESS","RESCAN_CA","RESCAN_CA_ITERATOR"]

    fig1 = plt.figure()
    boxplot1 = df.boxplot(column=all_columns ,rot=30, fontsize=10, figsize=(6,8))
    boxplot1.set_ylabel("COMPUTATION THROUGHPUT IN MB/s")
    boxplot1.set_xlabel("TESTED FUNCTIONS")
    fig1.tight_layout()
    fig1.savefig(f"{outfile}.png", format='png', dpi=300)

    fig2 = plt.figure()
    boxplot2 = df.boxplot(column=syncmer_columns, rot=30, fontsize=10, figsize=(6,8))
    boxplot2.set_ylabel("COMPUTATION THROUGHPUT IN MB/s")
    boxplot2.set_xlabel("TESTED FUNCTIONS")
    fig2.tight_layout()
    fig2.savefig(f"{outfile}_focus.png", format='png', dpi=300)


