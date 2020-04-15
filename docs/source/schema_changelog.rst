.. _dbchangelog:

MySQL database schema changelog
===============================

Below is a history of schema changes.
For maps of each schema version: :ref:`schema maps <schemamaps>`.


Schema version 9
****************

Substantial restructuring of trna and tmrna tables.

**Column moved**

    - trna_structures.Structure to trna.Structure

**Column removed**

    - trna.Sequence
    - trna.Product
    - trna.InfernalScore

**Column renamed**

    - trna.TrnaID TO trna.GeneID
    - tmrna.TmrnaID TO tmrna.GeneID

**Column modified**

    - trna.Structure datatype
    - trna.AminoAcid datatype
    - trna.PhageID position in table and foreign key constraint
    - trna.LocusTag position in table
    - tmrna.PhageID position in table and foreign key constraint
    - tmrna.LocusTag position in table

**Column added**

    - trna.Source
    - trna.Name
    - tmrna.Length
    - tmrna.Name

**Table removed**

    - trna_structures



Schema version 8
****************

Modified foreign key constraint.

**Column modified**

    - gene.PhamID foreign key constraint



Schema version 7
****************

Renamed several columns. Dropped the pham table and renamed the pham_color table.

**Column removed**

    - phage.Cluster
    - pham_color.ID (primary key)
    - gene.ID

**Column renamed**

    - phage.HostStrain TO phage.HostGenus
    - phage.Cluster2 TO phage.Cluster
    - phage.Subcluster2 TO phage.Subcluster
    - phage.SequenceLength TO phage.Length
    - pham_color.Name TO pham_color.PhamID



**Column modified**

    - pham_color.PhamID now primary key


**Column added**

    - gene.Parts
    - gene.PhamID


**Table removed**

    - pham


**Table renamed**

    - pham_color TO pham



Schema version 6
****************

Standardized column nomenclature.

**Column renamed**

    - domain.id TO domain.ID
    - domain.description TO domain.Description
    - gene_domain.id TO gene_domain.ID
    - gene_domain.expect TO gene_domain.Expect
    - phage.status TO phage.Status
    - pham.name TO pham.Name
    - pham_color.id TO pham_color.ID
    - pham_color.name TO pham_color.Name
    - pham_color.color TO pham_color.Color
    - version.version TO version.Version
    - gene.translation TO gene.Translation
    - gene.id TO gene.ID
    - gene.cdd_status TO gene.DomainStatus
    - version.schema_version TO version.SchemaVersion
    - domain.hit_id TO domain.HitID
    - gene_domain.hit_id TO gene_domain.HitID
    - gene_domain.query_start TO gene_domain.QueryStart
    - gene_domain.query_end TO gene_domain.QueryEnd



Schema version 5
****************

Removed several tables and columns.

**Table removed**

    - node
    - host_range
    - host
    - pham_history
    - pham_old
    - scores_summary

**Column removed**

    - phage.Prophage
    - phage.Isolated
    - phage.ProphageOffset
    - phage.DateLastSearched
    - phage.AnnotationQC
    - gene.StartCodon
    - gene.StopCodon
    - gene.GC1
    - gene.GC2
    - gene.GC3
    - gene.GC
    - gene.LeftNeighbor
    - gene.RightNeighbor
    - gene.clustalw_status
    - gene.blast_status
    - gene.TypeID
    - pham.orderAdded


**Column modified**

    - phage.status datatype = enum('unknown','draft','final')



Schema version 4
****************

Added several tables.

**Table created**

    - tmrna
    - trna
    - trna_structures

**Column modified**

    - gene.translation datatype = VARCHAR(5000)



Schema version 3
****************

Added/removed several columns.

**Column created**

    - gene.LocusTag
    - version.schema_version
    - phage.Subcluster2
    - phage.Cluster2


**Column removed**

    - phage.Program



Schema version 2
****************

Added several columns.

**Column created**

    - phage.AnnotationAuthor
    - phage.Program
    - phage.AnnotationQC
    - phage.RetrieveRecord



Schema version 1
****************

Misc changes to maintain referential integrity.


**Table created**

    - version

**Column created**

    - gene.cdd_status

**CASCADE setting updated**

    - gene.PhageID
    - gene_domain.GeneID
    - pham.GeneID
    - scores_summary.query
    - scores_summary.subject



Schema version 0
****************

The base schema.
