#!/bin/bash

##TODO:
#1. Add checkpointing

module load kentutils/0.0 lastz/1.04.15
full_path=$(realpath $0)
basepath=$(dirname ${full_path})

usage() { echo -e "\nUsage: $0 -h help -s <speciestable> -p <projectdir>\n\n" 1>&2; exit 1; }

no_args="true"
while getopts ":hs:p:" option; do
    case "${option}" in
				h) usage;;
        s) speciestable=${OPTARG};;
        p) projectdir=${OPTARG};;
	:) printf "missing argument for -%s\n" "$OPTARG" >&2; usage;;
	 \?) printf "illegal option: -%s\n" "$OPTARG" >&2; usage;;
        *) usage;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && usage
[[ -f "${speciestable}" ]] || usage
if [[ ! -d "${projectdir}" ]]
then
	mkdir -p "${projectdir}"
fi

while IFS= read -r line
do
	[[ "$line" =~ ^#.*$ ]] && continue
	IFS=$'\t' read -r -a sinfo <<< "$line"
	
	mkdir -p ${projectdir}/genomes/${sinfo[1]}

	if [[ ${sinfo[3]} == http* ]]
	then
		wget --continue -O ${projectdir}/genomes/${sinfo[1]}/$(basename ${sinfo[3]}) ${sinfo[3]}
	fi
	if [[ ${sinfo[3]} == ftp* ]]
	then
		wget --continue -O ${projectdir}/genomes/${sinfo[1]}/$(basename ${sinfo[3]}) ${sinfo[3]}
	fi
	if [[ ${sinfo[3]} == LOCAL:* ]]
	then
		fname=$(echo ${sinfo[3]} | sed 's/LOCAL://')
		cp $fname ${projectdir}/genomes/${sinfo[1]}/$(basename ${sinfo[3]})
	fi

	if [[ ${sinfo[3]} == *.gz ]]
	then
		fa=${projectdir}/genomes/${sinfo[1]}/$(basename ${sinfo[3]} .gz)
		gunzip -c $fa.gz >$fa
	else
		fa=${projectdir}/genomes/${sinfo[1]}/$(basename ${sinfo[3]})
	fi
	##convert to twobit format
	faToTwoBit $fa $fa.2bit
	##write lastz capsule file
	lastz_32 $fa[multiple] --writecapsule=$fa.capsule --ambiguous=iupac
	##create sizes file
	faSize -detailed $fa >$fa.sizes
	##split genome into 1Mb sequence chunks with five chunks in a file
	splitdir=$(dirname ${fa})/splits
	mkdir -p ${splitdir}
	perl ${basepath}/splitFasta.pl -i $fa -o $splitdir/`basename $fa .fna`.C -w 1000000 -l 1000000
done < "${speciestable}"
