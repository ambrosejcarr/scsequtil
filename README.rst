
Tools for easy parsing of Fastq and GTF files in python
=======================================================

Installation:
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
#. Printing pandas tables to pdf: ``scsequtil.table.print_pdf``
#. Plotting multiple identical figures, one per axis: ``scsequtil.plot.grid.AxesGrid``
#. Plotting scatterplots with categorical (``scsequtil.plot.scatter.categorical``) or
   continuous (``scsequtil.plot.scatter.continuous``) feature coloring

Under Development
-----------------
#. Package to automate spawning of sharedmemory enabled process pools
#. Package to automate spawning of multiple processes to execute a bash command
   containing pipes.