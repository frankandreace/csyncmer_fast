all: fig_throughput.pdf

parameter_sweep.tsv: run_parameter_sweep.sh
	bash run_parameter_sweep.sh > $@.tmp && mv $@.tmp $@ || { rm -f $@.tmp; false; }

data.tsv: run_benchmarks.sh
	bash run_benchmarks.sh > $@.tmp && mv $@.tmp $@ || { rm -f $@.tmp; false; }

fig_throughput.pdf: plot_throughput.py data.tsv syng.tsv
	python3 plot_throughput.py data.tsv -s syng.tsv -o fig_throughput.pdf

clean:
	rm -f data.tsv parameter_sweep.tsv fig_throughput.pdf fig_throughput.png

.PHONY: all clean
