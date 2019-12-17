.. _flatfile:

GenBank-formatted flat files
============================

GenBank-formatted flat files are a common data format used to represent an entire phage genome, and a detailed description of this data structure can be found at NCBI website :gbfff:`GenBank-formatted flat file <>`. This is a structured text file that systematically stores diverse types of information about the genome.

A flat file can be generated for any genome at any annotation stage using:

    1. GenBank
    2. :dnamaster:`DNA Master <>`
    3. :pecaan:`PECAAN <>`
    4. :biopython:`Biopython <>` (:ref:`Cock et al., 2009 <cock2009>`)

Flat file fields, such as LOCUS, DEFINITION, and REFERENCE-AUTHORS provide information regarding the entire record, while others, such as FEATURES, provide information about particular regions of the sequence in the record, such as tRNA or CDS genes. Data from flat files are stored in the *phage* and *gene* tables.

.. _figflatfile:

.. figure:: /images/flatfile_parsing.jpg

    Summary of how a GenBank-formatted flat file is parsed
