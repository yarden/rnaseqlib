.. include:: <isogrk3.txt>

.. documentation master file, created by .. 
   sphinx-quickstart on Fri Oct 22 16:50:57 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



.. contents::


``rnaseqlib``
=============

``rnaseqlib`` is a simple, lightweight pipeline for RNA-Seq analysis. 


Features
========

* Pipeline for analyzing RNA sequencing data:
  - Maps reads to genome and splice junctions
  - Computes basic quality control statistics
  - Outputs RPKM values for genes
* Supports:
  - mRNA sequencing (mRNA-Seq)
  - Ribosome profiling data (Ribo-Seq)
  - *CLIP for RNA-binding proteins (CLIP-Seq) [In progress]*
  - *SELEX-Seq [In progress]*


Updates
=======

**2012**

* ...


Installation
============

To install with `distribute`'s `easy_install`: ::

  easy_install rnaseqlib

For local installation, use: ::

  easy_install --user -U rnaseqlib

If you don't have `easy_install`, install it as an ordinary Python package: ::

  python setup.py install 

Or: ::

  python setup.py install --prefix=/your/local/dir

For local installation.


Running ``rnaseqlib``
=====================
.. _Running rnaseqlib

Initializing an RNA base for your genome
----------------------------------------

The first step is to create a set of files, called an *RNA Base*, required to run `rnaseqlib` for a particular genome. For example, ::

  rna_pipeline.py --init mm9 --output-dir pipeline_init

This will create a directory called `mm9` in `pipeline_init` containing the necessary files for running the pipeline on the ``mm9`` genome.
The initialization procedure will, among other things, do the following:

  * Download the genome sequence files from UCSC (in FASTA format)
  * Download gene tables from UCSC for Ensembl genes, UCSC knownGenes, and RefSeq genes.
  * Computate the coordinates for exons, introns, constitutive exons, constitutive exons in coding regions, and other 
    useful features of gene tables
  * Index the genome files using ``bowtie-build``


Configuration: specifying your data and mapping parameters
----------------------------------------------------------
.. _config:

The settings of the RNA-Seq pipeline are specified through a single settings file that contains four main sections:

  * ``pipeline``: what data type is used (e.g. RNA-Seq, Ribo-Seq, etc.)
  * ``pipeline-files``: where the pipeline initialization files are stored
  * ``mapping``: parameters related to mapping of the data (e.g. what mapper to use, where the genome index is.)
  * ``data``: where the input sequence files are and where the results should be outputted

The following is an example settings file for a set of mRNA-Seq samples: ::

  ##
  ## Pipeline settings for RNA-Seq
  ##
  [pipeline]
  # Data type to use (e.g. 'rnaseq' or 'riboseq')
  data_type = rnaseq

  [pipeline-files]
  init_dir = ./mm9

  [mapping]
  mapper = tophat
  bowtie_path = bowtie
  tophat_path = tophat
  tophat_index = mm9_index
  cluster_type = bsub
  genome = mm9
  tophat_options = --bowtie1 --min-anchor-length 4 
  tophat_gtf = ./mm9/ucsc/knownGene.gtf
  readlen = 40
  overhanglen = 4
  paired = True
  paired_end_frag = 300
  stranded = fr-first

  [data]
  indir = ./input_dir
  outdir = ./results

  sequence_files = [
      ["sample1_p1.fastq.gz", "sample1_p1"],
      ["sample1_p2.fastq.gz", "sample1_p2"],
      ["sample2_p1.fastq.gz", "sample2_p1"],
      ["sample2_p2.fastq.gz", "sample2_p2"]]

  sample_groups = [["sample1", ["sample1_p1", "sample1_p2"]],
                   ["sample2", ["sample2_p1", "sample2_p2"]]]

Running the pipeline
--------------------

To run the pipeline, use the ``--run`` option: ::

  rna_pipeline.py --run --settings ./settings.txt --output-dir ./my_results

where ``settings.txt`` is the pipeline settings file and ``my_results`` is a directory where the pipeline output should go.

 
Command-line options
--------------------

``run_pipeline.py`` is the main driver script of the pipeline. It takes the following main arguments: ::


  --run                 Run pipeline.

  --settings      
                        Settings filename.

  --init
                        Initialize the pipeline. Takes as input a genome, e.g.
                        mm9 or hg18.

  --output-dir
                        Output directory.



Frequently Asked Questions (FAQ)
================================

.. note:: 
  Section under construction

.. _refs:
.. _katz:


Authors
=======

``rnaseqlib`` was written by Yarden Katz.

Acknowledgements
================

Thanks to:

* Noah Spies (MIT)
* Eric Wang (MIT)



.. _MISO: http://genes.mit.edu/burgelab/miso/
.. _samtools: http://samtools.sourceforge.net/
.. _pysam: http://code.google.com/p/pysam/
.. _SAM: http://samtools.sourceforge.net/SAM1.pdf
.. _GFF: http://www.sequenceontology.org/gff3.shtml
.. _Drosophila melanogaster alternative events (modENCODE): http://genes.mit.edu/burgelab/miso/annotations/modENCODE_alt_events.zip
.. _Mouse genome (mm9) alternative events: http://genes.mit.edu/burgelab/miso/annotations/mm9_alt_events.zip
.. _Human genome (hg18) alternative events: http://genes.mit.edu/burgelab/miso/annotations/hg18_alt_events.zip
.. _Human genome (hg19) alternative events: http://genes.mit.edu/burgelab/miso/annotations/hg19_alt_events.zip
.. _Indexed mm9 annotations: http://genes.mit.edu/burgelab/miso/annotations/mm9/pickled/
.. _Indexed hg18 annotations: http://genes.mit.edu/burgelab/miso/annotations/hg18/pickled/
.. _Bowtie: http://bowtie-bio.sourceforge.net/
.. _Tophat: http://tophat.cbcb.umd.edu/
.. _IGV: http://www.broadinstitute.org/igv/
.. _PolyA DB: http://polya.umdnj.edu/polyadb/
.. _repository: https://github.com/yarden/MISO
.. _Perl script: http://seqanswers.com/forums/showthread.php?t=3201&highlight=GFF3
.. _miso-users: http://mailman.mit.edu/mailman/listinfo/miso-users
.. _White House adopts MISO: http://www.mediabistro.com/fishbowldc/white-house-soup-of-the-day-64_b53593#.Tp2c76k31tA.gmail

.. Indices and tables
.. ^^^^^^^^^^^^^^^^^^

.. * :ref:`genindex`
.. * :ref:`search`

