#!/bin/bash

basepath=`dirname $0`

input="$basepath/metadata/species.txt"
set -ex

while IFS= read -r line
do
	[[ "$line" =~ ^#.*$ ]] && continue
	IFS=$'\t' read -r -a sinfo <<< "$line"
	[[ ${sinfo[7]} -le 18 ]] && continue
	mkdir -p $basepath/../genomes/${sinfo[1]}
	wget --continue -O $basepath/../genomes/${sinfo[1]}/`basename ${sinfo[3]}` ${sinfo[3]}
	wget --continue -O $basepath/../genomes/${sinfo[1]}/`basename ${sinfo[6]}` ${sinfo[6]}
	fa=$basepath/../genomes/${sinfo[1]}/`basename ${sinfo[3]} .gz`
	gunzip -c $fa.gz >$fa
	##convert to twobit format
	$basepath/../software/kentutils/faToTwoBit $fa $fa.2bit
	##write lastz capsule file
	$basepath/../software/lastz-1.04.03/src/lastz_32 $fa[multiple] --writecapsule=$fa.capsule --ambiguous=iupac
	##create sizes file
	$basepath/../software/kentutils/faSize -detailed $fa >$fa.sizes
	##split genome into 1Mb sequence chunks with five chunks in a file
	splitdir=`dirname $fa`/splits
	mkdir -p $splitdir
	perl $basepath/splitFasta.pl -i $fa -o $splitdir/`basename $fa .fna`.C -w 1000000 -l 5000000
done < "$input"
