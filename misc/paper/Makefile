all: fig_throughput.pdf

parameter_sweep.tsv: run_parameter_sweep.sh
	bash run_parameter_sweep.sh > $@.tmp && mv $@.tmp $@ || { rm -f $@.tmp; false; }

data.tsv: run_benchmarks.sh
	bash run_benchmarks.sh > $@.tmp && mv $@.tmp $@ || { rm -f $@.tmp; false; }

# syng_syncmer.tsv (from syng_benchmark.sh) is optional: it splits panel (c) bars
SYNG_SYNCMER = $(wildcard syng_syncmer.tsv)

fig_throughput.pdf: plot_throughput.py data.tsv syng.tsv $(SYNG_SYNCMER)
	python3 plot_throughput.py data.tsv -s syng.tsv $(if $(SYNG_SYNCMER),--syng-syncmer $(SYNG_SYNCMER)) -o fig_throughput.pdf

clean:
	rm -f data.tsv parameter_sweep.tsv fig_throughput.pdf fig_throughput.png

.PHONY: all clean
