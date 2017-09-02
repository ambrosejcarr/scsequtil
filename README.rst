
Simple tools to ease parsing and construction of fastq and gtf files in python
==============================================================================

installation:
-------------
``pip3 install scsequtil``

Fastq:
------
``scsequtil.fastq.Reader``
``scsequtil.fastq.Record``

Gtf:
----
``scsequtil.gtf.Reader``
``scsequtil.gtf.Record``


Convenience Functions:
----------------------
#. printing pandas tables to pdf: ``scsequtil.table.print_pdf``
#. plotting multiple identical figures, one per axis: ``scsequtil.plot.grid.AxesGrid``
#. plotting scatterplots with categorical (``scsequtil.plot.scatter.categorical``) or
   continuous (``scsequtil.plot.scatter.continuous``) feature coloring

Under Development
-----------------
#. package to automate spawning of sharedmemory enabled process pools
#. package to automate spawning of multiple processes to execute a bash command
   containing pipes.