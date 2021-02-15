#!/bin/bash

full_path=$(realpath $0)
basepath=`dirname $full_path`

lastzbin="/g/data/te53/hrp561/wga/software/lastz-1.04.03/src/lastz_32"
input="$basepath/metadata/species.txt"
#set -ex

targetcapsules=()
splitsdir=()
while IFS= read -r line
do
	[[ "$line" =~ ^#.*$ ]] && continue
	IFS=$'\t' read -r -a sinfo <<< "$line"
	fa=$basepath/../genomes/${sinfo[1]}/`basename ${sinfo[3]} .gz`
	fa=$(realpath $fa)
	targetcapsules+=($fa.capsule)
	splitsdir+=(`dirname $fa`/splits)
	speciescodes+=(${sinfo[1]})
done < "$input"

for i in "${!targetcapsules[@]}"
do
	#echo ${targetcapsules[i]}
	for j in "${!targetcapsules[@]}"
	do
		[[ $j -le $i ]] && continue
		[[ $j -ne 18 ]] && continue
		alnoutput="${splitsdir[j]}/../aln/vs${speciescodes[i]}"
		mkdir -p $alnoutput
		alnoutput=$(realpath $alnoutput)
#		echo $alnoutput
		for s in ${splitsdir[j]}/*.fa
		do
			outputpath=$alnoutput/`basename $s .fa`.vs.${speciescodes[i]}
			outputbase=`basename $s .fa`.vs.${speciescodes[i]}
			echo "sh $basepath/runcmd.sh \"$lastzbin --targetcapsule=${targetcapsules[i]} $s K=2400 L=3000 Y=9400 H=2000 --ambiguous=iupac --format=axt --output=$outputpath.axt\" $outputbase.lastz $outputpath.lastz.done 0"
		done
	done
done
