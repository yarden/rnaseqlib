.. include:: <isogrk3.txt>

.. documentation master file, created by .. 
   sphinx-quickstart on Fri Oct 22 16:50:57 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. toctree::

.. contents::
  :depth: 3


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

  - CLIP for RNA-binding proteins (CLIP-Seq) 

  - SELEX-Seq (Bind-n-Seq)

* Contains utilities for processing `MISO`_ output


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


Updates
=======

**2012**

* Released ``rnaseqlib``


Installation
============

Get ``rnaseqlib`` from the GitHub repository (http://github.com/yarden/rnaseqlib). To install with `distribute`_'s ``easy_install``: ::

  easy_install rnaseqlib

For local installation with `distribute`_, use: ::

  easy_install --user -U rnaseqlib

If you don't have ``easy_install``, install it as an ordinary Python package: ::

  python setup.py install 

Or for local installation: ::

  python setup.py install --prefix=/your/local/dir

Dependencies
------------

Requires `Bedtools`_ and `Samtools`_, and several commonly used Python modules
(like `pandas`_, all of which are installed automatically
using the Python package manager.) The initialization step requires Jim Kent's 
`genePredToGtf`_ utility for processing UCSC genome browser tables.



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
  * Compute the coordinates for exons, introns, constitutive exons, constitutive exons in coding regions, and other 
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

Example MISO pipeline
---------------------

The following pipeline steps start with a settings file describing a set of samples and run MISO on the samples, generating summary files and comparisons for relevant samples, and combining the results in unfiltered and filtered (based on read counts) form.

1. **Run** MISO on samples: ::

  misowrap.py run misowrap_settings.txt miso_output/

The raw MISO output will be placed in the directory defined in ``misowrap_settings.txt`` and the logs of the run will be placed in ``miso_output/`` directory.

2. Once run is completed, **summarize samples**: ::

  misowrap.py summarize misowrap_settings.txt miso_output/

This summarizes samples and places summary files for each sample in relevant place as inferred by sample organization in ``misowrap_settings.txt``. Logs are placed in ``miso_output/`` directory.

3. **Compare samples** to each other. Uses the ``sample_groups`` parameter in settings file to determine which 
pairwise comparisons to run. ::

  misowrap.py compare misowrap_settings.txt miso_output/

This compare samples to each other and places the output of comparisons in a ``comparisons`` subdirectory of the MISO output directory given in ``misowrap_settings.txt``. Logs are placed in the ``miso_output/`` directory.

4. **Filter comparisons** based on reads counts, using the count filters defined for each event type in ``misowrap_settings.txt``: ::

  misowrap.py filter-comparisons misowrap_settings.txt miso_output/

This creates filtered versions of the comparison ``.miso_bf`` files and places them in the directory ``filtered_events`` in the MISO comparisons directory (which is the ``comparisons`` directory in the MISO output directory given in ``misowrap_settings.txt``). For example: ::

  - MISO output directory (defined by settings file)
    - comparisons/
        - filtered_events/
            - SE/
              - SE samples here...
              - ...
            - RI/ 
              - RI samples here...
              - ...

5. **Combine comparisons** into single files and place them in the directory where our comparisons are. Takes the ``comparisons`` directory (which contains MISO comparisons) and the settings file ::

  misowrap.py combine-comparisons miso_output/comparisons/ misowrap_settings.txt --output-dir miso_output/comparisons/

This outputs single files that contain information about events pooled from all the samples described in the settings file. It also creates a filtered version of these combined events where the read filters (defined for each event) in ``misowrap_settings.txt`` are applied to the samples. The resulting output structure is: ::

  - MISO output directory (defined by settings file)
    - comparisons/
      - combined_comparisons/    # combined comparisons (unfiltered)
      - filtered_events/         # filtered events
        - combined_comparisons/  # combined comparisons (filtered)
          - SE.miso_bf           # combined, filtered comparisons for all SE events
          - RI.miso_bf           # combined, filtered comparisons for all RI events
          - ...


Creating custom GFF annotations
===============================

``rnaseqlib`` has a set of scripts in the ``gff`` module that can generate a GFF annotation of exon-centric alternative events which can be used by MISO for quantitation using RNA-Seq. Given a set of transcripts in the UCSC genePred format, these scripts will build a "splice graph" representation of how each exon in the transcript can be spliced, and then produce a set of possible events (categorized into alternative event types) by traversing this graph. For example, suppose we have a genePred table called ``ensGene.txt``: ::

  $ head -n 3 ensGene.txt
  585     ENST00000619216.1       chr1    -       17368   17436   17368   17368   1       17368,  17436,  0       MIR6859-2       none    none    -1,
  585     ENST00000473358.1       chr1    +       29553   31097   29553   29553   3       29553,30563,30975,      30039,30667,31097,      0       MIR1302-11      none    none    -1,-1,-1,
  585     ENST00000469289.1       chr1    +       30266   31109   30266   30266   2       30266,30975,    30667,31109,    0       MIR1302-11      none    none    -1,-1,

We can produce a set of exon-centric GFF events from this table using the ``gff_make_annotation`` script from ``rnaseqlib``: ::

  $ gff_make_annotation ./ ./gff --flanking-rule commonshortest --genome-label hg38
  Making GFF alternative events annotation...
    - UCSC tables read from: ./
    - Output dir: ./gff
  Loaded 1 UCSC tables.
  Loading tables...
  Populating graph...
  Adding splice edges from table ensGene.txt
  Populating graph took 37.50 seconds
  Reading table ./ensGene.txt
  Generating skipped exons (SE)
  Generating mutually exclusive exons (MXE)
  Generating alternative 3' splice sites (A3SS)
  Generating alternative 5' splice sites (A5SS)
  Outputting retained introns...
  Defining retained introns (RI)
  Took 1.60 minutes to make the annotation.

The call to ``gff_make_annotation`` says: load all genePred tables in the current directory (``./``) and output the set of GFF events into the directory ``./gff``. The ``--flanking-shortest`` parameter specifies how to pick the flanking exons of an alternative event when there are multiple options (e.g. which flanking exons to use for an alternatively skipped exon.) In the above call, we chose to the take the common shortest region as flanking exons. The ``--genome-label`` option specifies what the name of the GFF annotation should be. Note that ``rnaseqlib`` will only load genePred tables if they are named ``ensGene.txt``, ``refGene.txt`` or ``knownGene.txt``. Any genePred files with these names that are in the input directory will be parsed and used to populate the splice graph, so that information from multiple annotation genePred tables can be used to make the GFF events.

This call should generate an output directory with a GFF file for each of the event types: ::

  $ ls ./gff/commonshortest/
  A3SS.hg38.gff3  A5SS.hg38.gff3  MXE.hg38.gff3  RI.hg38.gff3  SE.hg38.gff3


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

