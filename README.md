# annotation-database-constructor

This program takes as input files from the `gene-check` pipeline and produces a `sqlite3` database. This database is constructed so that each column represents one **feature** found in the data. Each row represents one alignment, and is uniquely keyed on the alignment ID number. The database contains a table for each 1-to-1 comparison between a reference transcriptome and the target genome.

To run this program, you need to have cython installed. This is required for the `twobit` python library, which will be automatically compiled and installed by the makefile. This also means you need to have ability to add modules to your python installation. If you want to do this manually, look in `lib/twobit/` for the pyx file.

Once this is setup, run the following commands:

```git submodule update --init
make all```

Then to run the program, you need to modify the makefile to point to the right files, then run either `make run` to run the annotation pipeline. You can also use `make details` to build the secondary database from the gene-check-details files.

The features currently represented in this are:

1. geneID - the Gencode ID for the gene this transcript is from.
2. Gene - the common name for this transcript.
3. sourceChrom - the chromosome the gene is on in the source genome.
4. sourceStart - the 0-based start location of this transcript in the source genome.
5. sourceEnd - the 0-based stop location of this transcript in the source genome.
6. sourceStrand - the strand the transcript is on in the source genome.
7. destChrom - the chromosome this alignment is on in the target genome.
8. destStart -the 0-based start location of this alignment in the target genome.
9. destEnd -the 0-based stop location of this alignment in the target genome.
10. sourceStrand - the strand the alignment is on in the target genome.
11. beginStart - what position in the transcript does the CDS start? -1 if the CDS does not start with an ATG.
12. endStop - are the last three bases of the CDS in this alignment a stop codon? 1 if True, 0 if False.
13. inFrameStop - is there a in-frame stop codon in the CDS of this alignment? 1 if True, 0 if False.
14. noCds - does this alignment lack a CDS? Defined as having thickStart/thickStop. 1 if True, 0 if False.
15. cdsGap - does there exist a intron between CDS that is too short? Too short is currently defined as <=30bp. 1 if True, 0 if False.
16. cdsMult3Gap - does there exist a intron between CDS that is both too short and a multiple of 3. 1 if True, 0 if False.
17. utrGap - does there exist a intron between non-coding exons that is too short? Too short is currently defined as <= 30bp. 1 if True, 0 if False.
18. cdsUnknownSplice - does there exist a intron beween CDS whose splice sites do not fit one of the known sites, `GT..AG`, `GC..AG`, `AT..AC`? 1 if True, 0 if False.
19. cdsNonCanonSplice?  does there exist a intron beween CDS whose splice sites do not fit the canonical splice site, `GT..AG`? 1 if True, 0 if False.
20. utrUnknownSplice - does there exist a intron beween non-coding exons whose splice sites do not fit one of the known sites, `GT..AG`, `GC..AG`, `AT..AC`? 1 if True, 0 if False.
21. utrNonCanonSplice - does there exist a intron beween non-coding exons whose splice sites do not fit the canonical splice site, `GT..AG`? 1 if True, 0 if False.
22. alignmentCoverage - a boolean value between 0 and 1 that represents the alignment coverage, calculated as (matches + mismatches) / (matches + mismatches + query_insertions). 
23. alignmentIdentity - a boolean value between 0 and 1 that represents alignment identity, calculated as matches / (matches + mismatches + query_insertions).
24. alignmentAbutsLeft - does this alignment hit or overlap with the left edge of a assembly scaffold in the target genome? 1 if True, 0 if False.
25. alignmentAbutsRight - does this alignment hit or overlap with the right edge of a assembly scaffold in the target genome? 1 if True, 0 if False.
26. unknownBases - a integer value representing the number of Ns (unknown bases) in this alignment.
27. numberScaffoldGap - a integer value representing the number of scaffold gaps (represented as 100 Ns) in this alignment.
28. alignmentPartialMap - Does the query transcript map entirely? 1 if True, 0 if False.
