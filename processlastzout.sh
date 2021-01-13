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
		alnoutput=$(realpath "${splitsdir[j]}/../aln/vs${speciescodes[i]}")
		pairbase=$alndir/${speciescodes[j]}.vs.${speciescodes[i]}
		mergedout=$pairbase.axt
		mergedout_h=$pairbase.headers
		mergedout_a=$pairbase.alignments
		chkpnt=$pairbase.lastz.chkpnt

		rm -f $mergedout $mergedout_a $mergedout_h $chkpnt

		for s in ${splitsdir[j]}/*.fa
		do
			outputpath=$alnoutput/`basename $s .fa`.vs.${speciescodes[i]}
			outputbase=`basename $s .fa`.vs.${speciescodes[i]}
			if [[ -e "$outputpath.lastz.done" && "`tail -n1 $outputpath.lastz.done | cut -f3 -d ','`" == " EXIT_STATUS:0" ]]
			then
				perl -lne 'print $_ if ($_ =~ /^#/)' $outputpath.axt >> $pairbase.headers
				perl -lne 'print $_ unless ($_ =~ /^#/)' $outputpath.axt >> $pairbase.alignments
				tail -n4 $outputpath.lastz.done >> $pairbase.lastz.chkpnt
			else
				echo "Either $outputpath.lastz.done does not exist or task had non-zero exit status"
			fi
		done
		##merge headers and alignments
		mv $pairbase.headers $pairbase.lastz.axt
		perl -lne 'if ($_=~/^\d+/) { @a = split (" ", $_); $a[4] =~ /(\S+)\.s(\d+)\.e\d+/; $c=$1;$s=$2; $a[4] = $c; $a[5] = $s + $a[5]; $a[6] = $s + $a[6]; $a[0] = $counter ? $counter : 0; $counter++; print join (" ", @a); } else { print $_ }' $pairbase.alignments >> $pairbase.lastz.axt
		rm -f $pairbase.alignments
		echo sh $basepath/runcmd.sh \"sh $basepath/getMAF.sh $pairbase ${twobits[i]} ${twobits[j]} ${sizes[i]} ${sizes[j]} ${speciescodes[i]} ${speciescodes[j]}\" ${speciescodes[j]}.vs.${speciescodes[i]}.chainnet $pairbase.chainnet.done 0 >>$mafcmd
		#rm -rf $alnoutput
	done
done

