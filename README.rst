
# Simple tools to ease parsing and construction of fastq and gtf files in python

# installation:
pip3 install scsequtil

# fastq:
scsequtil.fastq.Reader
scsequtil.fastq.Record

# gtf
scsequtil.gtf.Reader
scsequtil.gtf.Record

# also contains some convenience functions for

1. printing pandas tables to pdf: scsequtil.table.print_pdf
2. plotting multiple identical figures, one per axis: scsequtil.plot.grid.AxesGrid
3. plotting scatterplots with categorical (scsequtil.plot.scatter.categorical) or
   continuous (scsequtil.plot.scatter.continuous) feature coloring

# under development
1. package to automate spawning of sharedmemory enabled process pools
2. package to automate spawning of multiple processes to execute a bash command
   containing pipes.