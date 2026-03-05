#!/usr/bin/env python3
"""Generate publication-quality syncmer throughput bar charts."""

import argparse
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Method display names and colors (order = bar order left→right)
# Slot 0 reserved for future "digest" bar
METHODS = [
    ("digest",          "digest"),
    ("seqhash",         "seqhash (syng)"),
    ("rescan",          "Rescan (C)"),
    ("simd-minimizers", "simd-minimizers (Rust)"),
    ("twostack",        "csyncmer_fast (single)"),
    ("multi-8",         "csyncmer_fast (multi-8)"),
]
METHOD_ORDER = {m[0]: i for i, m in enumerate(METHODS)}
COLORS  = ["#999999", "#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#e377c2"]  # grey, red, blue, green, orange, pink
ALPHAS  = [0.4, 0.4, 0.4, 0.4, 1.0, 1.0]  # csyncmer_fast methods at full opacity

# Methods to hide from the plot (data still loaded, just not drawn)
HIDE = {"seqhash", "rescan"}

DATASETS = [
    ("chm13", "(a) CHM13 (FASTA)"),
    ("hifi",  "(b) HiFi reads (FASTQ)"),
]

# syng panel: map config names in syng.tsv → METHODS keys
SYNG_MAP = {"seqhash": "seqhash", "csyncmer": "twostack", "multi8": "multi-8"}
SYNG_ORDER = ["seqhash", "csyncmer", "multi8"]  # bar order left→right


def read_data(path):
    """Read data.tsv → {dataset: {method: throughput_gbps}}."""
    data = defaultdict(dict)
    with open(path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            dataset, method, throughput = line.split("\t")
            data[dataset][method] = float(throughput)
    return data


def read_syng(path):
    """Read syng.tsv → {config: median_s}."""
    data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            data[parts[0]] = float(parts[1])
    return data


def plot(data, syng_data, outbase):
    FS = 8  # base font size — tune this single value to scale all text
    plt.rcParams.update({
        "font.size": FS,
        "axes.labelsize": FS,
        "axes.titlesize": FS + 1,
        "xtick.labelsize": FS - 0.5,
        "ytick.labelsize": FS - 0.5,
        "legend.fontsize": FS - 1,
    })

    visible = [(m, c, a) for (m, c, a) in zip(METHODS, COLORS, ALPHAS) if m[0] not in HIDE]
    n_visible = len(visible)
    bar_width = 0.36
    group_width = n_visible * bar_width
    group_gap = group_width * 0.45  # space between dataset groups

    has_syng = bool(syng_data)
    ncols = 2 if has_syng else 1
    fig, axes = plt.subplots(1, ncols, figsize=(4.5 if has_syng else 3.3, 1.8),
                             gridspec_kw={"width_ratios": [1.8, 0.65]} if has_syng else None)
    if ncols == 1:
        axes = [axes]
    fig.subplots_adjust(wspace=0.25)

    ax_tp = axes[0]  # throughput axis (panels a + b combined)

    max_y = 0
    for ds_key, _ in DATASETS:
        ds_data = data.get(ds_key, {})
        for m_key, _ in METHODS:
            if m_key in ds_data and m_key not in HIDE:
                max_y = max(max_y, ds_data[m_key])
    max_y *= 1.18

    # Plot both dataset groups on the same axis
    group_centers = []
    for g, (ds_key, ds_label) in enumerate(DATASETS):
        ds_data = data.get(ds_key, {})
        group_center = g * (group_width + group_gap)
        group_centers.append(group_center)
        for i, ((m_key, m_label), color, alpha) in enumerate(visible):
            v = ds_data.get(m_key)
            if v is None:
                continue
            x = group_center + i * bar_width - (n_visible - 1) * bar_width / 2
            ax_tp.bar(x, v, width=bar_width * 0.95, color=color, alpha=alpha,
                      edgecolor="black", linewidth=0.3, label=m_label if g == 0 else None)
            ax_tp.text(x, v + max_y * 0.01, f"{v:.2f}", ha="center",
                       va="bottom", fontsize=FS - 2)

    ax_tp.set_xticks(group_centers)
    ax_tp.set_xticklabels([dl for _, dl in DATASETS], fontsize=FS - 0.5)
    ax_tp.tick_params(axis="x", length=0, pad=3)
    ax_tp.spines["top"].set_visible(False)
    ax_tp.spines["right"].set_visible(False)
    ax_tp.tick_params(axis="y", length=2, pad=1)
    ax_tp.set_ylabel("Throughput (GB/s)", fontsize=FS, labelpad=4)
    if max_y > 0:
        ax_tp.set_ylim(0, max_y)
    # ── Panel (c): syng execution time ──
    if has_syng:
        ax_syng = axes[1]
        # Build lookup from method key → (color, alpha, label)
        method_style = {m[0]: (c, a, m[1]) for m, c, a in zip(METHODS, COLORS, ALPHAS)}
        n_syng = len(SYNG_ORDER)
        max_t = max(syng_data.get(c, 0) for c in SYNG_ORDER) * 1.18

        for i, cfg in enumerate(SYNG_ORDER):
            m_key = SYNG_MAP[cfg]
            color, alpha, m_label = method_style[m_key]
            v = syng_data.get(cfg)
            if v is None:
                continue
            x = i * bar_width - (n_syng - 1) * bar_width / 2
            ax_syng.bar(x, v, width=bar_width * 0.95, color=color, alpha=alpha,
                        edgecolor="black", linewidth=0.3, label=m_label)
            ax_syng.text(x, v + max_t * 0.01, f"{v:.1f}s", ha="center",
                         va="bottom", fontsize=FS - 2)

        ax_syng.set_xticks([])
        ax_syng.set_xlabel("(c) syng (8 threads)", fontsize=FS - 0.5, labelpad=2)
        ax_syng.spines["top"].set_visible(False)
        ax_syng.spines["left"].set_visible(False)
        ax_syng.yaxis.tick_right()
        ax_syng.yaxis.set_label_position("right")
        ax_syng.set_ylabel("Time (s)", fontsize=FS, labelpad=6)
        ax_syng.spines["right"].set_position(("outward", 5))
        ax_syng.tick_params(axis="y", length=2, pad=1)
        ax_syng.set_ylim(0, max_t)

    # Single legend above all panels — merge handles from all axes
    seen = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in seen:
                seen[l] = h
    handles = list(seen.values())
    labels = list(seen.keys())
    fig.legend(handles, labels, loc="upper center", ncol=3,
               frameon=False, fontsize=FS - 1, bbox_to_anchor=(0.5, 1.09),
               handlelength=1.0, handletextpad=0.3, columnspacing=0.8)

    # Save PDF + PNG
    pdf_path = outbase if outbase.endswith(".pdf") else outbase + ".pdf"
    png_path = pdf_path.replace(".pdf", ".png")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.12)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.12)
    print(f"Saved {pdf_path} and {png_path}", file=sys.stderr)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", help="Path to data.tsv")
    parser.add_argument("-s", "--syng", default=None,
                        help="Path to syng.tsv (adds panel c)")
    parser.add_argument("-o", "--output", default="fig_throughput.pdf",
                        help="Output file (default: fig_throughput.pdf)")
    args = parser.parse_args()

    data = read_data(args.data)
    syng_data = read_syng(args.syng) if args.syng else {}
    plot(data, syng_data, args.output)


if __name__ == "__main__":
    main()
