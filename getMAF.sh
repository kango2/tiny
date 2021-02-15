#!/bin/bash

set -e

##load Jim Kent's Utilities available at http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ with download instructions
#module load kentutils

pairbase=$1
target2bit=$2
query2bit=$3
targetsize=$4
querysize=$5
targetscode=$6
queryscode=$7
codebase=$8

##fix lastz alignments as coordinates for reverse strand were wrongly recorded
perl $codebase/fixlastz.pl $pairbase.lastz.axt.gz $pairbase.lastz.corrected.axt $querysize
##chain lastz alignments
axtChain -minScore=3000 -linearGap=medium $pairbase.lastz.corrected.axt $target2bit $query2bit $pairbase.chains &>/dev/null
##sort chains
chainSort $pairbase.chains $pairbase.sorted.chains
rm -f $pairbase.chains
##chainPreNet
chainPreNet $pairbase.sorted.chains $targetsize $querysize $pairbase.prenet.chains
##net chains
chainNet $pairbase.prenet.chains $targetsize $querysize $pairbase.target.tmpnet $pairbase.query.tmpnet
##target netSyntenic
netSyntenic $pairbase.target.tmpnet $pairbase.target.net
##query netSyntenic
netSyntenic $pairbase.query.tmpnet $pairbase.query.net
##netToaxt
#netToAxt $pairbase.target.net $pairbase.prenet.chains $target2bit $query2bit $pairbase.net.axt
##axtsort
#axtSort $pairbase.net.axt $pairbase.sorted.net.axt
##axtToMaf
#axtToMaf -tPrefix="$targetscode"- -qPrefix="$queryscode"- $pairbase.sorted.net.axt $targetsize $querysize $pairbase.maf
#zgrep -v ^# $i | perl -lne 'if ($_=~/net (\S+) (\d+)/) { $tc=$1;$tcl=$2 } elsif ($_=~/^ fill (\d+) (\d+) (\S+) (\S) (\d+) (\d+)/) { $ts=$1; $tl=$2; $qc=$3; $strand = $4; $qs=$5; $ql=$6; $te=$ts+$tl; $qe=$qs+$ql; print "$tc\t$ts\t$te\t$qc\t$qs\t$qe\t$strand" }' >/g/data/te53/hrp561/wga/lastzaln/`basename $i .net.gz`.chainpairs.tab
##remove unwanted files
rm -f $pairbase.net.axt $pairbase.query.tmpnet $pairbase.target.tmpnet $pairbase.sorted.chains $pairbase.prenet.chains
##gzip wanted files
gzip --force $pairbase.lastz.corrected.axt $pairbase.query.net $pairbase.target.net
