#!/bin/bash

set -e

export PATH=/g/data/te53/hrp561/wga/software/kentutils/:$PATH

pairbase=$1
target2bit=$2
query2bit=$3
targetsize=$4
querysize=$5
targetscode=$6
queryscode=$7

axtChain -minScore=3000 -linearGap=medium $pairbase.lastz.axt $target2bit $query2bit $pairbase.chains &>/dev/null
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
netToAxt $pairbase.target.net $pairbase.prenet.chains $target2bit $query2bit $pairbase.net.axt
##axtsort
axtSort $pairbase.net.axt $pairbase.sorted.net.axt
##axtToMaf
axtToMaf -tPrefix="$targetscode"- -qPrefix="$queryscode"- $pairbase.sorted.net.axt $targetsize $querysize $pairbase.maf
rm -f $pairbase.net.axt $pairbase.query.tmpnet $pairbase.target.tmpnet 

