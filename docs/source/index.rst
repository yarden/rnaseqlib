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

* Supports several sequencing experiments:

  - mRNA sequencing (mRNA-Seq)

  - Ribosome profiling data (Ribo-Seq)

  - CLIP for RNA-binding proteins (CLIP-Seq) *(In progress)*

  - SELEX-Seq *(In progress)*
* Contains utilities for processing `MISO`_ output


Updates
=======

**2012**

* Released ``rnaseqlib``


Installation
============

Get ``rnaseqlib`` from the GitHub repository (http://github.com/yarden/rnaseqlib). To install with `distribute`_'s ``easy_install``: ::

  easy_install rnaseqlib

For local installation, use: ::

  easy_install --user -U rnaseqlib

If you don't have `easy_install`, install it as an ordinary Python package: ::

  python setup.py install 

Or for local installation with distribute: ::

  python setup.py install --prefix=/your/local/dir

Dependencies
------------

Requires `Bedtools`_ and `Samtools`_, and several commonly used Python modules
(like `pandas`_, all of which are installed automatically
using the Python package manager.) The initialization step requires Jim Kent's 
`genePredToGtf`_ utility for processing UCSC genome browser tables.


Design principles
=================

``rnaseqlib`` follows three simple design principles. It is intended to be:

  1. **Simple:** provides minimalistic support for RNA-Seq. Performs only simple computations that 
     are applicable to nearly all experiments -- complexities that are specific to certain
     experiments/libraries are left as post-processing steps for the user.


  2. **Lightweight:** minimal dependencies. Relies mostly on Python and commonly used
     genomic packages (such as Bedtools), to avoid software bloat and complex installation.


  3. **Compact:** produces and consumes compressed files, so that it can be used in projects 
     with hundreds of samples.


Running ``rnaseqlib``
=====================

There are two steps to running ``rnaseqlib``: First, creating a set of initialization
files for your genome (called an *RNA Base*) -- this is done once per genome.
Second, writing a configuration file that describes your samples and library parameters,
which can then be used to run the pipeline.

Initializing an RNA base for your genome
----------------------------------------

The first step is to create a set of files, called an *RNA Base*, 
required to run ``rnaseqlib`` for a particular genome. For example, ::

  rna_pipeline.py --init mm9 --output-dir pipeline_init

This will create a directory called ``mm9`` in ``pipeline_init`` containing the 
necessary files for running the pipeline on the mm9 genome.
The initialization procedure will, among other things, do the following:

  * Download the genome sequence files from UCSC (in FASTA format)
  * Download gene tables from UCSC for Ensembl genes, UCSC knownGenes, and RefSeq genes.
  * Computate the coordinates for exons, introns, constitutive exons, constitutive exons in coding regions, and other 
    useful features of gene tables
  * Index the genome files using ``bowtie-build``

For the mouse and human genomes, sequences of ribosomal RNA (rRNA) are automatically
downloaded from NCBI and built into the Bowtie index as ``chrRibo``. The pipeline
relies on ``chrRibo`` in later steps to filter out rRNA reads and measure the level
of rRNA contamination in libraries.


Configuration: specifying your data and mapping parameters
----------------------------------------------------------
.. _config:

The settings of the RNA-Seq pipeline are specified through a single settings file 
that contains four main sections:

  * ``pipeline``: what data type is used (e.g. RNA-Seq, Ribo-Seq, etc.)
  * ``pipeline-files``: where the pipeline initialization files are stored
  * ``mapping``: parameters related to mapping of the data (e.g. what mapper to use, 
    where the genome index is.)
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

where ``settings.txt`` is the pipeline settings file and ``my_results`` is a 
directory where the pipeline output should go.

An example settings file and small FASTQ files for an mRNA-Seq dataset are available 
in the ``examples/rnaseq`` directory of the pipeline. The paths in the settings 
file ``examples/rnaseq/rnaseq_settings.txt`` need to be edited to reflect the paths of 
various files on your own filesystem (e.g. ``init_dir`` needs to be set to the 
location of your RNA base.) The rawdata files for the examples are two paired-end 
mRNA-Seq samples are available in the ``examples/rnaseq/fastq`` directory.
 
Command-line options
--------------------

``run_pipeline.py`` is the main driver script of the pipeline. 
It takes the following main arguments: ::


  --run                 Run pipeline.

  --settings      
                        Settings filename.

  --init
                        Initialize the pipeline. Takes as input a genome, e.g.
                        mm9 or hg18.

  --output-dir
                        Output directory.




Settings file features and options
----------------------------------

Below are the main parameters to be used in the ``rnaseqlib`` settings file. 
Parameters that have default values are considered optional and can be omitted 
from the file.

* ``mapper``: The mapper to use (currently, only ``tophat`` is supported.)

* ``bowtie_path``: Path to ``bowtie`` program (optional). Default is ``"bowtie"``, meaning ``bowtie`` must be on your path.

* ``tophat_path``: Path to ``tophat`` program (optional). Default is ``"tophat"``, meaning ``tophat`` must be on your path.

* ``tophat_index``: Path to genome index that should be used by ``tophat`` mapping. This index is typically built by ``bowtie-build`` and is automatically made as part of the RNA Base using the ``--init`` option. 

* ``cluster_type``: The ``cluster_type`` settings under the ``[mapping]`` section of the settings file specifies how distinct samples should be processed. It can be set to one of three values:

  * ``bsub``: Run distinct samples as separate jobs on a cluster using ``bsub``
  
  * ``qsub``: Same as ``bsub``, but using the ``qsub`` submission system

  * ``none``: Run the pipeline using multiple cores on the local computer (does not make use of a cluster.)

* ``genome``: Genome to be used, e.g. ``mm9`` or ``hg18``

* ``tophat_options``: String of command-line options to be passed to ``tophat`` when it is invoked for mapping. E.g. ``tophat_options = --bowtie1 --min-anchor-length 4``, to signal to Tophat to use ``bowtie1`` for mapping and use a minimum of 4 base overhang for junction reads. This string is appended to the Tophat call and so must contain valid Tophat arguments for the call to succeed.

* ``tophat_gtf``: GTF file of known gene models to be used with Tophat (optional). Used to tell Tophat about known junctions.

* ``readlen``: Read length of input BAM files.

* ``overhanglen``: Overhang restriction on junctions (optional). Default is 1.

* ``paired``: Whether the data is paired-end or not. If paired-end, specify ``True``, if single-end, specify ``False``.

* ``paired_end_frag``: Average fragment length (optional). Only used for paired-end runs. Used internally as an argument to Tophat to specify the expected fragment length.

* ``stranded``: If data set is strand-specific, specify the strand convention (optional). Uses the same strand conventions as Tophat (e.g. ``fr-first``).

Creating and processing MISO output with ``misowrap``
=====================================================

.. note::

  Section under construction

``misowrap`` is a utility for running `MISO`_ and processing its output for a set of samples.
It takes a configuration file as described above that lists a set of BAM files and their sample 
labels, and automatically runs MISO on these samples and generates pairwise comparisons between them.

It can also create a set of filtered events, that are selected to meet read coverage criteria.
Filtered events are outputted to a table containing various gene features from the UCSC tables
generated by the pipeline's RNA base.


Features
--------

``misowrap`` takes the following arguments: ::


  --run                 Run MISO on a set of events. Takes a settings
                        filename.
  --summarize
                        Run MISO summarize on a set of samples. Takes a
                        settings filename.

  --compare             Run MISO sample comparisons on all pairwise
                        comparisons. Takes a settings filename.

  --filter              Filter a set of MISO events. Takes a settings
                        filename.

  --compute-insert-lens
                        Compute insert lengths for a set of BAM files. takes a
                        settings filename.

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
.. _bedtools: http://code.google.com/p/bedtools/
.. _genePredToGtf: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _pandas: http://pandas.pydata.org/
.. _distribute: http://packages.python.org/distribute/
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

