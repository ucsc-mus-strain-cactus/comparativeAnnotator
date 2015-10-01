TODO: THIS IS OUT OF DATE

# comparative annotation pipeline

This pipeline is used to create comparative annotation of aligned whole genomes. It does this by taking a set of target genomes and one reference genome whose annotations will be mapped over, then checked for correctness. An assemblyHub is produced automatically to visualize the results.


# INSTALLATION
To install this program, you need the following things.
1. sqlite3 in your path with version 3.8.7.4 or above
2. The python package `pyfaidx` which can be gotten through pip


`annotation-database-constructor` constructs `sqlite3` databases from alignments and BED/genePred files from the `transmap` and `gene-check` pipeline. Thus, the input for each genome is:

1. PSL of alignments where the query is a transcript and the target is the genome the annotations are being transferred to.
2. BED file representing the `transmap` output where the automatically annotated new transcript is.

In addition, an attributes file mapping the unique transcript ID names to attributes such as common gene name is used.

Note that these PSLs/BED files should be uniquely keyed, meaning that alignments should have a unique number added to it.

Given these inputs, and 2bit files representing source and target *genome*, this pipeline constructs three complementary databases:

1. **classify** - this database has values of 1 for True and 0 for False in all cells. Represents boolean classifications of each transcript alignment for the categories below.
2. **details** - this database has the same columns as classify, but with details represented as a string of BED records.
3. **attributes** - this database stores attributes about the transcripts such as their gene name.

If you want to add more classifiers, add them to the `classifiers.py` script in src/.


1. CodingInsertions - are there insertions that are not a multiple of 3 in coding sequence?
2. CodingMult3Insertions - same as CodingInsertions, but only multiples of 3.
3. CodingDeletions - are there deletions that are not a multiple of 3 in coding sequence?
4. CodingMult3Deletions - same as CodingDeletions, but only multiples of 3.
5. Rearrangements - looks for jumps in PSL coordinates that are indicative of rearrangements. This is defined as a indel happening and then the coordinates going the other direction.
6. FrameMismatch - Frameshifts are caused by coding indels that are not a multiple of 3.
7. AlignmentAbutsLeft - does this alignment hit or overlap with the left edge of a assembly scaffold in the target genome?
8. AlignmentAbutsRight - does this alignment hit or overlap with the right edge of a assembly scaffold in the target genome?
9. AlignmentPartialMap - Does the query transcript NOT map entirely?
10. BadFrame - is the CDS a multiple of 3?
11. EndStop - does the CDS start with 'ATG'?
12. CdsGap - does there exist an intron between CDS exons that is too short? Too short is currently defined as <=30bp. Only reports such gaps if they are not a multiple of 3.
13. CdsMult3Gap - same as CdsGap but reports only multiples of 3.
14. UtrGap - same as CdsGap, but for UTR introns.
15. CdsUnknownSplice - does there exist a intron beween CDS whose splice sites do not fit one of the known sites, `GT..AG`, `GC..AG`, `AT..AC`?
19. CdsNonCanonSplice - does there exist a intron beween CDS whose splice sites do not fit the canonical splice site, `GT..AG`?
16. UtrUnknownSplice - does there exist a intron beween non-coding exons whose splice sites do not fit one of the known sites, `GT..AG`, `GC..AG`, `AT..AC`?
17. UtrNonCanonSplice - does there exist a intron beween non-coding exons whose splice sites do not fit the canonical splice site, `GT..AG`?
18. EndStop - does the CDS end with a stop codon? ('TAA', 'TGA', 'TAG')
19. InFrameStop - is there a stop codon within the coding frame?
20. NoCds - is there no annotated CDS?
21. ScaffoldGap - Does this alignment span a scaffold gap (represented as 100 Ns)?
22. UnknownBases - Are there Ns in the alignment?
23. UnknownCdsBases - same as UnknownBases, but only if the Ns are in the CDS is this true.
24. Nonsynonymous - looks for nonsynonymous mutations. Does not report mutations in frameshifted regions.
25. Synonymous - looks for synonymous mutations. Does not report mutations in frameshifted regions.
26. Paralogy - Does this query transcript have more than one target alignment?
