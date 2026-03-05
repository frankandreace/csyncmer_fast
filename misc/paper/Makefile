all: fig_throughput.pdf

data.tsv: run_benchmarks.sh
	bash run_benchmarks.sh > data.tsv

fig_throughput.pdf: plot_throughput.py data.tsv
	python3 plot_throughput.py data.tsv -o fig_throughput.pdf

clean:
	rm -f data.tsv fig_throughput.pdf fig_throughput.png

.PHONY: all clean
