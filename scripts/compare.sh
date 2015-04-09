#!/bin/bash

grep 'UtrGap' $2 | awk '{print $1, $2, $3}' OFS="\t" > my.bed
grep 'utrGap' $1 > mark.bed
echo 'track name="UtrGap_mySubMark"' > $3/UtrGap.my_sub_mark.bed
echo 'track name="UtrGap_markSubMy"' > $3/UtrGap.mark_sub_my.bed
subtractBed -a my.bed -b mark.bed >> $3/UtrGap.my_sub_mark.bed
subtractBed -b my.bed -a mark.bed >> $3/UtrGap.mark_sub_my.bed

grep 'CdsGap\|CdsMult3Gap' $2 | awk '{print $1, $2, $3}' OFS="\t" > my.bed
grep 'cdsGap' $1 > mark.bed
echo 'track name="CdsGap_mySubMark"' > $3/CdsGap.my_sub_mark.bed
echo 'track name="CdsGap_markSubMy"' > $3/CdsGap.mark_sub_my.bed
subtractBed -a my.bed -b mark.bed >> $3/CdsGap.my_sub_mark.bed
subtractBed -b my.bed -a mark.bed >> $3/CdsGap.mark_sub_my.bed

grep 'BeginStart' $2 | awk '{print $1, $2, $3}' OFS="\t" > my.bed
grep 'noStart' $1 > mark.bed
echo 'track name="BeginStart_mySubMark"' > $3/BeginStart.my_sub_mark.bed
echo 'track name="BeginStart_markSubMy"' > $3/BeginStart.mark_sub_my.bed
subtractBed -a my.bed -b mark.bed >> $3/BeginStart.my_sub_mark.bed
subtractBed -b my.bed -a mark.bed >> $3/BeginStart.mark_sub_my.bed

grep 'EndStop' $2 | awk '{print $1, $2, $3}' OFS="\t" > my.bed
grep 'noStop' $1 > mark.bed
echo 'track name="NoStop_mySubMark"' > $3/EndStop.my_sub_mark.bed
echo 'track name="NoStop_markSubMy"' > $3/EndStop.mark_sub_my.bed
subtractBed -a my.bed -b mark.bed >> $3/EndStop.my_sub_mark.bed
subtractBed -b my.bed -a mark.bed >> $3/EndStop.mark_sub_my.bed


grep 'CdsUnknownSplice' $2 | awk '{print $1, $2, $3}' OFS="\t" > my.bed
grep 'unknownCdsSplice' $1 > mark.bed
echo 'track name="CdsUnKnownSplice_mySubMark"' > $3/CdsUnKnownSplice.my_sub_mark.bed
echo 'track name="CdsUnKnownSplice_markSubMy"' > $3/CdsUnKnownSplice.mark_sub_my.bed
subtractBed -a my.bed -b mark.bed >> $3/CdsUnKnownSplice.my_sub_mark.bed
subtractBed -b my.bed -a mark.bed >> $3/CdsUnKnownSplice.mark_sub_my.bed

rm my.bed mark.bed