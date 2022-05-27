.. _phamerate:

phamerate: create gene phamilies
================================

Phameration is the process whereby phamilies ("phams") are constructed from groups of putatively homologous
protein-coding genes. This procedure allows us to distill the genetic diversity in a set of genomes down to its most
discrete units, which in turn facilitates examination of the evolutionary relationships between very distantly related
genomes.

The ``pdm_utils phamerate`` tool provides access to a pham assembly pipeline using MMseqs2
(:ref:`Steinegger and Söding, 2017 <bibliography>`), a homology detection suite purpose built for searching and
clustering large genomics datasets.

To phamerate gene products using the MMseqs2 pipeline::

    > python3 -m pdm_utils phamerate Actino_Draft

**The pipeline will crash if the mmseqs2 software suite has not been installed or is not globally executable**

The argument 'Actino_Draft' indicates the name of the database whose gene products are to be phamerated.

This pipeline has been optimized from the ground up to construct phamilies of globally homologous phage proteins. In
general we expect that functions can be propagated across members within the same pham.

Pham colors (for display in Phamerator) are inserted into the *pham* table of the database. Any phams that are
unchanged (or now include one or more newly added genes) between rounds of phameration will have their pham
designation and color preserved.

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