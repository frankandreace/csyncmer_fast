#!/usr/bin/env python3

import pandas as pd
import matplotlib as mtpl
from sys import argv



if __name__ == "__main__" :
    
    file = argv[1]
    outfile = argv[2]

    df = pd.read_csv(file,delimiter="\t")

    # boxplot = df.boxplot(column=["HASHING","NAIVE","DEQUE","RESCAN","RESCAN_CIRCULAR_ARRAY","RESCAN_CA_ITERATOR"] ,rot=30, fontsize=10, figsize=(8,8))
    # mtpl.pyplot.savefig(f"{outfile}.png")

    boxplot = df.boxplot(column=["NAIVE","DEQUE","RESCAN","RESCAN_CIRCULAR_ARRAY","RESCAN_CA_ITERATOR"] ,rot=30, fontsize=10, figsize=(8,8))
    mtpl.pyplot.savefig(f"{outfile}_focus.png", format='png')


