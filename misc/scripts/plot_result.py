#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from sys import argv


if __name__ == "__main__" :
    
    file = argv[1]
    outfile = argv[2]

    df = pd.read_csv(file,delimiter="\t")

    all_columns = ["SYNGH_JUST_HASHING","NTH_JUST_HASHING","SYNGH_RESCAN_NO_ITERATOR","SYNGH_RESCAN_ITERATOR","NTH_ITERATOR","SYNGH_DURBIN_ITERATOR","SYNGH_DEQUE","NTH_DEQUE","SYNGH_NAIVE","SYNGH_RESCAN_CA","SYNGH_RESCAN_CA_ITERATOR","SYNGH_RESCAN_CA_BRANCHLESS"]

    reversed_col = list(reversed(all_columns))

    missing_columns = [col for col in all_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"The following columns are not in the DataFrame: {missing_columns}")


    syncmer_columns = ["SYNGH_DURBIN_ITERATOR","SYNGH_RESCAN_ITERATOR","SYNGH_RESCAN_NO_ITERATOR","SYNGH_JUST_HASHING"]

    my_figsize=(8, 6)

    fig1 = plt.figure(figsize=my_figsize)
    boxplot1 = df.boxplot(vert=False, column=reversed_col ,rot=0, fontsize=8)
    boxplot1.set_xlabel("COMPUTATION THROUGHPUT IN MB/s")
    boxplot1.set_ylabel("TESTED FUNCTIONS")
    boxplot1.set_title("Speed Benchmark. NTH = NT Hash; SYNGH = Syng Hash.")
    fig1.tight_layout()
    fig1.savefig(f"{outfile}.png", format='png', dpi=300)

    fig2 = plt.figure(figsize=my_figsize)
    boxplot2 = df.boxplot(vert=False, column=syncmer_columns, rot=0, fontsize=8)
    boxplot2.set_xlabel("COMPUTATION THROUGHPUT IN MB/s")
    boxplot2.set_ylabel("TESTED FUNCTIONS")
    boxplot2.set_title("Speed Benchmark. NTH = NT Hash; SYNGH = Syng Hash.")
    fig2.tight_layout()
    fig2.savefig(f"{outfile}_focus.png", format='png', dpi=300)


