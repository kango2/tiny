#!/bin/bash

set -ex

full_path=$(realpath $0)
basepath=`dirname $full_path`

lastzbin="/g/data/te53/hrp561/wga/software/lastz-1.04.03/src/lastz_32"
export PATH=$PATH:/g/data/te53/hrp561/wga/software/kentutils

input="$basepath/metadata/species.txt"
alndir=$(realpath "$basepath/../lastzaln")
mkdir -p $alndir

mafcmd=$alndir/getmaf.cmds.txt
rm -rf $mafcmd

targetcapsules=()
splitsdir=()
twobits=()
sizes=()

while IFS= read -r line
do
	[[ "$line" =~ ^#.*$ ]] && continue
	IFS=$'\t' read -r -a sinfo <<< "$line"
	fa=$basepath/../genomes/${sinfo[1]}/`basename ${sinfo[3]} .gz`
	fa=$(realpath $fa)
	targetcapsules+=($fa.capsule)
	twobits+=($fa.2bit)
	sizes+=($fa.sizes)
	splitsdir+=(`dirname $fa`/splits)
	speciescodes+=(${sinfo[1]})
done < "$input"

for i in "${!targetcapsules[@]}"
do
	for j in "${!targetcapsules[@]}"
	do
		[[ $j -le $i ]] && continue
		pairbase=$alndir/${speciescodes[j]}.vs.${speciescodes[i]}
		echo sh $basepath/runcmd.sh \"sh $basepath/getMAF.sh $pairbase ${twobits[i]} ${twobits[j]} ${sizes[i]} ${sizes[j]} ${speciescodes[i]} ${speciescodes[j]}\" ${speciescodes[j]}.vs.${speciescodes[i]}.chainnet $pairbase.chainnet.done 0 >>$mafcmd
	done
done

