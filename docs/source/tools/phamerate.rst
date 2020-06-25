.. _phamerate:

phamerate: create gene phamilies
================================

Phameration is the process whereby phamilies ("phams") are constructed from groups of putatively homologous
protein-coding genes. This procedure allows us to distill the genetic diversity in a set of genomes down to its most
discrete units, which in turn facilitates examination of the evolutionary relationships between very distantly related
genomes.

The ``pdm_utils phamerate`` tool provides access to two pipelines for phameration. The first pipeline uses MMseqs2
(:ref:`Steinegger and Söding, 2017 <bibliography>`), a homology detection suite purpose built for searching and
clustering large genomics datasets. The second pipeline mimics a pipeline that seems to be popular in the phage
genomics world, which uses blastp (:ref:`Altschul et al, 1990 <bibliography>`) for homology detection and Markov
Clustering (:ref:`van Dongen & Abreu-Goodger, 2012 <bibliography>`) to interpret the blastp output into clusters.

To phamerate gene products using the MMseqs2 pipeline::

    > python3 -m pdm_utils phamerate mmseqs Actinobacteriophage

**The pipeline will crash if the mmseqs software suite has not been installed or is not globally executable**

To phamerate gene products using the blastp --> MCL pipeline::

    > python3 -m pdm_utils phamerate blast-mcl Actinobacteriophage

**The pipeline will crash if the NCBI's blast+ package OR Markov Clustering have not been installed or are not
globally executable**

The argument 'Actinobacteriophage' indicates the name of the database whose gene products are to be phamerated.

Both pipelines have optional arguments that can be used to adjust the thresholds/algorithms used for clustering. The
mmseqs pipeline has been optimized from the ground up to construct phamilies of globally homologous phage proteins. In
general we expect that functions can be propagated across members within the same pham. The blast-mcl pipeline has not
been optimized nearly as rigorously, but it uses parameters I've seen commonly published for similar pipelines used
for phage genomics.

Regardless of which pipeline is used, new phams and their colors (for display in Phamerator) are inserted into the
*pham* table of the database. Any phams that are unchanged (or now include one or more newly added genes) between
rounds of phameration will have their pham designation and color preserved.

Notes for mmseqs pipeline
*************************

The mmseqs pipeline will run in two steps by default:

1.  Sequence-sequence clustering to construct pre-phamily HMM profiles
2.  Profile-consensus clustering to merge pre-phamilies into more sensitive phamilies

The second step can be skipped by using the --skip-hmm argument at the command line. Doing so without also adjusting
the parameters used in the first step will result in phams with relatively weak sensitivity, but few (if any) false
positives.

As previously mentioned, this pipeline's default parameters have been optimized to construct phams of globally
homologous phage proteins. Aside from minor tweaks to these parameters, in general the only reason to deviate from
the defaults would be to construct phams for a different use cases than comparing genomes on the basis of shared
genes or propagating functions. For example, if one's goal is to examine intragenic mosaicism, the default coverage
threshold is too high to identify most domain-linked sequences. In this case it's probably simpler to use the
blast-mcl pipeline with a low (or no) coverage cutoff.