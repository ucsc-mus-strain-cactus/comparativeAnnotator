# comparativeAnnotator

This pipeline attempts to classify and diagnose transcripts mapped over to target genomes using transMap. These transcripts are then cleaned up using a special form of Augustus, and finally a consensus gene set is produced.

See the main pipeline repo https://github.com/ucsc-mus-strain-cactus/pipeline


# DEPENDENCIES

1. sqlite3 with >= 3.8.7.4
2. python with the following packages: `pyfaidx`, `matplotlib`, `numpy`, `pandas`
3. R with the package `pvclust` (if you want to cluster classifiers)

For the full pipeline, you will also need the full Kent code base.


# INFORMATION

`comparativeAnnotator` constructs `sqlite3` databases from transMap results, attempting to diagnose issues with the transcripts. See the pipeline repo for the makefile wrappers that drive this program.

This pipeline produces a total of 3 databases in transMap mode, and an additional 3 in Augustus mode:

1. **classify** - This database contains the results of the classifiers described below.
2. **details** - This database has the same columns as classify, but with details represented as a string of BED records. These records will be loaded into the final assemblyHub.
3. **attributes** - this database stores attributes about the transcripts such as their gene name as well as some alignment metrics.

See the file etc/config.py for some example queries against these databases that are used to both determine the quality of a transcript as well as attempt to bin transcripts as likely having problems related to alignment or assembly.

# CLASSIFIERS
Each classifier is binary unless stated otherwise.

Reference classifiers. These classifiers do not depend on an alignment and so can be run against any single gene set + genome.

1. StartOutOfFrame - is the first coding base out of frame, as defined by the genePred exonFrames field?
2. BadFrame - is the lifted over CDS both longer than 25 amino acids and not a multiple of 3?
3. BeginStart - are the first 3 CDS bases 'ATG'?
4. EndStop - are the last 3 CDS bases one of 'TAA', 'TGA', 'TAG'?
5. CdsGap - are any of the CDS introns shorter than or equal to 30bp and are not a multiple of 3? Reports the number of such gaps.
6. CdsMult3Gap - are any of the CDS introns shorter than or equal to 30bp and are a multiple of 3? Reports the number of such gaps.
7. UtrGap - are any of the non-CDS introns shorter than or equal to 30bp? Reports the number of such gaps.
8. UnknownGap - are any of the introns shorter than or equal to 30bp, and do these introns contain Ns? Reports the number of such gaps.
9. CdsNonCanonSplice - are any of the CDS introns longer than 30bp and not of the canonical form 'GT:AG'? Reports the number of such introns.
10. CdsUnknownSplice - are any of the CDS introns longer than 30bp and not of any of the canonical or non-canonical splice sites ('AG:GC', 'GC:AG', 'AT:AC')? Reports the number of such introns.
11. CdsNonCanonSplice - are any of the non-CDS introns longer than 30bp and not of the canonical form 'GT:AG'? Reports the number of such introns.
12. CdsUnknownSplice - are any of the non-CDS introns longer than 30bp and not of any of the canonical or non-canonical splice sites ('AG:GC', 'GC:AG', 'AT:AC')? Reports the number of such introns.
13. SpliceContainsUnknownBases - are any of the introns longer than 30bp and have a unknown base in either the donor or acceptor sequence? Reports the number of such introns.
14. InFrameStop - in the frame the transcript is supposed to have (taking into account exonFrames), are there any in frame stop codons? Reports the number of in frame stops.
15. ShortCds - is the lifted over CDS shorter than or equal to 25 amino acids?
16. UnknownBases - how many bases in the lifted over mRNA sequence are Ns?
17. UnknownCdsBases - how many bases in the lifted over CDS are Ns?

Alignment classifiers. These classifiers rely on having a PSL mapping between a source and target transcript, as produced by transMap.
18. AlnExtendsOffContig - does the alignment go off the end of a contig? This is defined as the source transcript having unaligned portions while the alignment hits the edge of the contig.
19. AlnAbutsUnknownBases - does the 5' or 3' end of the transcript touch Ns?
20. HasOriginalIntrons - how many original intron boundaries (in transcript space) are not present in this lifted over transcript? This classifier looks at where the source transcript had splice sites in mRNA coordinates and compares this to the splice sites in the lifted over transcript. Reports the number of missing introns.
21. CodingInsertions - does the alignment introduce insertions to the target genome that are not a multiple of 3? Only evaluates transcripts with at least 25 amino acids. Reports the number of insertions.
22. CodingMult3Insertions - does the alignment introduce insertions to the target genome that are a multiple of 3? Only evalutes transcripts with at least 25 amino acids. Reports the number of insertions.
23. CodingDeletions - does the alignment introduce deletions to the target genome that are not a multiple of 3? Only evalutes transcripts with at least 25 amino acids. Reports the number of deletions.
24. CodingMult3Deletions - does the alignment introduce deletions to the target genome that are a multiple of 3? Only evalutes transcripts with at least 25 amino acids. Reports the number of deletions.
25. FrameShift - are there non multiple of 3 indels that lead to frameshifts? In details mode, attempts to determine if frame is restored by a subsequent indel. Reports the number of frame shifts.
26. AlignmentPartialMap - does the source transcript not map over entirely?
27. HasOriginalStart - based on the alignment, is the lifted over start codon aligned to the source start codon? This will be false when an incomplete alignment extends past the UTR to the CDS.
28. HasOriginalStop - based on the alignment, is the lifted over stop codon aligned to the source stop codon?
29. Nonsynonymous - are there any nonsynonymous changes in this transcript? Reports the number of nonsynonymous changes.
30. Synonymous - are there any synonymous changes in this transcript? Reports the number of synoynmous changes.
31. Paralogy - did transMap map this transcript to more than one place? Reports the number of extra alignments (n-1).

Augustus classifiers. These classifiers apply to transcripts resulting from AugustusTMR.
32. AugustusNotSameStrand - did augustusTMR produce a transcript on the same strand as the source transcript?
33. AugustusParalogy - did augustusTMR produce exactly one transcript for this source transcript?
34. AugustusExonGain - did augustusTMR add new exons relative to the source transcript?
35. AugustusExonLoss - did augustusTMR remove any exons?
36. AugustusNotSimilarInternalExonBoundaries - did augustus move any internal exons more than 30bp, after merging any transMap alignment gaps shorter than 50bp?
37. AugustusNotSimilarTerminalExonBoundaries - did augustus move any terminal exon boundaries more than 200bp, after merging any transMap alignment gaps shorter than 100bp?
38. AugustusNotSameStart - does the augustusTMR transcript have the exact same CDS start coordinates?
39. AugustusNotSameStop - does the augustusTMR transcript have the exact same CDS stop coordinates?
