#!/bin/bash

export PATH=/g/data/te53/hrp561/wga/software/kentutils/:$PATH

cd /g/data/te53/hrp561/wga/genomes/

for i in *_assembly_report.txt; 
do
	echo processing $i
	qbase=`dirname $i`/`basename $i _assembly_report.txt`
	qfbase=`basename $i _assembly_report.txt`
	tfbase=GCF_000002315.6_GRCg6a
	scode=`grep $i speciescode.txt | cut -f1`
	seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f7 | grep -v na`;
	if [[ ${#seqids} -eq 0 ]]; 
	then 
		seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f5 | grep -v na`;
	fi; 
	cmd="chainMergeSort"; 
	for j in $seqids; 
	do 
		if [ -e chrseq/$j.chain ]; 
		then 
			cmd="$cmd chrseq/$j.chain"; 
		fi;
	done;
	##merge and sort chromosomal chains
	eval $cmd >$qbase.chains;
	##chainPreNet
	chainPreNet $qbase.chains $tfbase.sizes $qbase.sizes $qbase.chains.prenet;
	##net chains
	chainNet $qbase.chains.prenet $tfbase.sizes $qbase.sizes $tfbase.vs.$qfbase.ttmpnet $qbase.qtmpnet
	##target netSyntenic
	netSyntenic $tfbase.vs.$qfbase.ttmpnet $tfbase.vs.$qfbase.tnet
	##query netSyntenic
	netSyntenic $qbase.qtmpnet $qbase.qnet
	##netToaxt
	netToAxt $tfbase.vs.$qfbase.tnet $qbase.chains.prenet "$tfbase"_genomic.2bit "$qbase"_genomic.2bit $qbase.axt
	##axtsort
	axtSort $qbase.axt $qbase.axtsorted
	##axtToMaf
	axtToMaf -tPrefix=CHICK- -qPrefix="$scode"- $qbase.axtsorted $tfbase.sizes $qbase.sizes $qbase.vs.CHICK.maf
done
