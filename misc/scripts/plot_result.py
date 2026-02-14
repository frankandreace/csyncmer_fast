#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from sys import argv


if __name__ == "__main__":

    if len(argv) < 3:
        print(f"Usage: {argv[0]} <benchmark.tsv> <output_prefix>")
        exit(1)

    file = argv[1]
    outfile = argv[2]

    df = pd.read_csv(file, delimiter="\t")

    # Get all columns from the TSV
    all_columns = list(df.columns)

    # Filter by prefix
    hashing_columns = [col for col in all_columns if col.startswith("HASH_")]
    syncmer_columns = [col for col in all_columns if col.startswith("SYNCMER_")]

    my_figsize = (10, max(6, len(all_columns) * 0.4))

    # Plot all implementations
    fig1 = plt.figure(figsize=my_figsize)
    boxplot1 = df.boxplot(vert=False, column=list(reversed(all_columns)), rot=0, fontsize=8)
    boxplot1.set_xlabel("Throughput (MB/s)")
    boxplot1.set_ylabel("Implementation")
    boxplot1.set_title("Speed Benchmark - All Implementations")
    fig1.tight_layout()
    fig1.savefig(f"{outfile}.png", format='png', dpi=300)
    print(f"Saved: {outfile}.png")

    # Plot syncmer implementations only
    if syncmer_columns:
        fig2 = plt.figure(figsize=(10, max(6, len(syncmer_columns) * 0.4)))
        boxplot2 = df.boxplot(vert=False, column=list(reversed(syncmer_columns)), rot=0, fontsize=8)
        boxplot2.set_xlabel("Throughput (MB/s)")
        boxplot2.set_ylabel("Implementation")
        boxplot2.set_title("Speed Benchmark - Syncmer Implementations")
        fig2.tight_layout()
        fig2.savefig(f"{outfile}_syncmers.png", format='png', dpi=300)
        print(f"Saved: {outfile}_syncmers.png")
