
# Micro-chromosome conservation in vertebrates

## Data preparation

Process was automated using `preparegenomes.sh` script. Basic steps were as follows.

1. Genome assemblies were downloaded from relevant sources (NCBI or DNAZoo). Check `metadata/species.txt` for details of species and download paths.
2. `.2bit`, `.capsule`, `.sizes` files were created.
3. Genome files were split into smaller chunks (1Mb sequence size and ~5Mb of total sequence per file) without overlaps. **This step allows for embarassingly parallel lastz alignments to capitalise on large HPC facilities.**

## One-way all-vs-all whole genome alignments

**Credits**

1. Workflow ideas for lastz alignments were borrowed from [Daren Card](https://github.com/darencard) (thanks mate) available [here](https://darencard.net/blog/2019-11-01-whole-genome-alignment-tutorial/).
2. Lastz alignment parameters were obtained from several sources:  
    * The bird genome alignment paper GitHub repo [here](https://github.com/gigascience/paper-zhang2014/blob/master/Whole_genome_alignment/pairwise/bin/lastz_CNM.pl).

    ``` 
    K=2400 L=3000 Y=9400 H=2000
    ```
    * Alignment sensitivity tests performed in [Sharma et al, 2017 paper](https://doi.org/10.1093/nar/gkx554).

    ```
    K = 2400, L = 3000, Y = 9400, H = 2000 for placental mammals
    K = 2400, L = 3000, Y = 3400, H = 2000 for non-placental mammals
    K = 1500, L = 2500 and W = 5  to find co-linear alignments in the un-aligning regions that are flanked by local alignments (gaps in the chains)
    ```
    * Ensembl Compara LastZ pairwise alignment settings for the GitHub Repo [here](https://github.com/Ensembl/ensembl-compara/blob/23bcb7ecaed4b6ea3251b22b1405d9d9e0d817bc/modules/Bio/EnsEMBL/Compara/PipeConfig/Lastz_conf.pm)

    ```
    default => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac', # ensembl genomes settings
    # Vertebrata
    7742    => 'T=1 K=3000 L=3000 H=2200 O=400 E=30 --ambiguous=iupac',
    # Catarrhini, Sus, Carnivora, Triticeae
    9526    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 Q='
    9822    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac',
    33554   => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac',
    147389  => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac --identity=75..100',
    # Vigna, Solanaceae
    3913    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    4070    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    #4107    => 'K=5000 L=5000 H=3000 O=400 E=30 --ambiguous=iupac M=10 --notransition --step=20',
    #4107    => 'K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac --notransition --step=20',
    ```
    * Default parameters of LastZ
    ```
    # hsp_threshold (K)      = 3000
    # gapped_threshold (L)   = 3000
    # x_drop (X)             = 910
    # y_drop (Y)             = 9400
    # gap_open_penalty (O)   = 400
    # gap_extend_penalty (E) = 30
    ```

*We used `K=2400 L=3000 Y=9400 H=2000 --ambiguous=iupac` as parameters for LastZ alignments for the project for all species.*

**List of species included in analyses**
|Species|Species Code|Common Name|NCBI Taxonomy ID|
|---|---|---|---|
|Gallus gallus|CHICK|Chicken|9031|
|Trachemys scripta elegans|TRASE|Red-eared slider turtle|31138|
|Gallus gallus|CHICKYO|Chicken Yeonsan Ogye|9031|
|Corvus monedula|CORMO|Jackdaw|30423|
|Cygnus olor|CYGOL|Mute swan|8869|
|Naja naja |NAJNA|Indian cobra|35670|
|Ornithorhynchus anatinus|ORNAN|Platypus|9258|
|Gopherus evgoodei|GOPEV|Goodes thornscrub tortoise|1825980|
|Dermochelys coriacea|DERCO|Leatherback sea turtle|27794|
|Chelonia mydas|CHEMY|Green sea turtle|8469|
|Lacerta agilis|LACAG|Sand lizard|80427|
|Thamnophis elegans|THAEL|Western terrestrial garter snake|35005|
|Podarcis muralis|PODMU|Common wall lizard|64176|
|Zootoca vivipara|ZOOVI|Common lizard|8524|
|Crotalus viridis|CROVV|Prairie rattlesnake|8742|
|Calypte anna|CALAN|Anna's hummingbird|9244|
|Falco rusticolus|FALRU|Gyrfalcon|120794|
|Aquila chrysaetos chrysaetos|AQUCH|Golden eagle|223781|
|Dromaius novaehollandiae|DRONO|Emu|8790|
|Salvator merianae|SAMER|Argentine black and white tegu|96440|
|Python bivittatus|PYTBI|Burmese python|176946|
|Gymnogyps californianus|GYMCA|California condor|33616|
|Alligator mississippiensis|ALLMI|American alligator|8496|
|Alligator sinensis|ALLSI|Chinese alligator|38654|
|Homo sapiens|HUMAN|Human|9606|
|Sarcophilus harrisii|SARHA|Tasmanian devil|9305|
|Branchiostoma floridae|BRAFL|Florida lancelet|7739|


```
##Download genomes
for i in `cut -f3 metadata.txt | grep -v fasta`; do wget --continue $i; done
##Download assembly report
for i in `cut -f6 metadata.txt | grep -v assemblyreport`; do wget --continue $i; done
##Create sequence sizes file
for i in *.fna.gz; do gunzip -c $i | perl -lne 'if ($_ =~ />(\S+)/){ $h=$1; push (@seq, $h); } else { $l{$h}+=length($_) } END { map { print "$_\t$l{$_}" } @seq }' >`dirname $i`/`basename $i _genomic.fna.gz`.sizes; done
##Unzip genomes
for i in *.fna.gz; do zcat $i > `basename $i .gz`; done
##index genome files
for i in *.fna; do samtools faidx $i; done
##get sequence names for the query species
for i in *_assembly_report.txt; do seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f7 | grep -v na`; if [[ ${#seqids} -eq 0 ]]; then seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f5 | grep -v na`; fi; for j in $seqids; do if [ ! -e chrseq/$j.fa ]; then samtools faidx `basename $i _assembly_report.txt`_genomic.fna $j >chrseq/$j.fa; faToTwoBit chrseq/$j.fa chrseq/$j.2bit; fi; done; done
```


Create target capsule file
```
/g/data/te53/hrp561/wga/software/lastz-1.04.03/src/lastz_32 /g/data/te53/hrp561/wga/genomes/GCF_000002315.6_GRCg6a_genomic.fna[multiple] --writecapsule=/g/data/te53/hrp561/wga/genomes/GCF_000002315.6_GRCg6a_genomic.capsule --ambiguous=iupac
```

Create lastz commands: takes time to run per chromosome
```
( for i in /g/data/te53/hrp561/wga/genomes/chrseq/*.fa; do echo sh /g/data/te53/phase2_20200312/utils/runcmd.sh \"/g/data/te53/hrp561/wga/software/lastz-1.04.03/src/lastz_32 --targetcapsule=/g/data/te53/hrp561/wga/genomes/GCF_000002315.6_GRCg6a_genomic.capsule $i K=2400 L=3000 Y=9400 H=2000 --ambiguous=iupac --format=axt --output=$i.axt\" `basename $i .fa`.lastz `dirname $i`/`basename $i .fa`.lastz.done 0; done ) >lastzcmds.txt
```

Launch commands
```
qsub -j oe -o pbslogs/ -l walltime=48:00:00,ncpus=48,mem=190GB -N laln -V -v commandsfile=lastzcmds.txt,ncpupercmd=1,joblog=pbslogs/lastz.parallel.0.log runcmdsparallel.sh
```
Things to organise
```
##axtChain, default score used is 5000 in the gitlab script,
##histogram shows the number of alignments increases sharply at 3000 score
##a graph to show this would be good per species.
##needs twoBit file
for i in genomes/chrseq/*.fa; do echo processing $i; software/kentutils/faToTwoBit $i `dirname $i`/`basename $i .fa`.2bit; done

##following runs very quick and it is to chain lastZ alignments
for i in /g/data/te53/hrp561/wga/genomes/chrseq/*.2bit; do d=`dirname $i`/`basename $i .2bit`.lastz.done;  if [ "`tail -n1 $d | cut -f3 -d','`" == " EXIT_STATUS:0" ]; then /g/data/te53/hrp561/wga/software/kentutils/axtChain -minScore=3000 -linearGap=medium `dirname $i`/`basename $i .2bit`.fa.axt /g/data/te53/hrp561/wga/genomes/GCF_000002315.6_GRCg6a_genomic.2bit $i `dirname $i`/`basename $i .2bit`.chain; fi; done


##merge and sort chains for each species
##Bird genome paper doesn't have the following steps which are listed in Daren's blog
##patchChain, RepeatFiller, chainCleaner

for i in *_assembly_report.txt; do seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f7 | grep -v na`; if [[ ${#seqids} -eq 0 ]]; then seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f5 | grep -v na`; fi; cmd="/g/data/te53/hrp561/wga/software/kentutils/chainMergeSort"; for j in $seqids; do if [ -e chrseq/$j.chain ]; then cmd="$cmd /g/data/te53/hrp561/wga/genomes/chrseq/$j.chain"; fi; done; eval $cmd >`dirname $i`/`basename $i _assembly_report.txt`.chains; /g/data/te53/hrp561/wga/software/kentutils/chainPreNet `dirname $i`/`basename $i _assembly_report.txt`.chains /g/data/te53/hrp561/wga/genomes/GCF_000002315.6_GRCg6a.sizes `dirname $i`/`basename $i _assembly_report.txt`.sizes `dirname $i`/`basename $i _assembly_report.txt`.chains.prenet; done

../software/kentutils/chainNet GCA_009769625.1_bCygOlo1.pri.prenet.chains GCF_000002315.6_GRCg6a.sizes GCA_009769625.1_bCygOlo1.pri.sizes GCF_000002315.6_GRCg6a.vs.GCA_009769625.1_bCygOlo1.pri.ttmpnet GCA_009769625.1_bCygOlo1.pri.qtmpnet

../software/kentutils/netSyntenic GCF_000002315.6_GRCg6a.vs.GCA_009769625.1_bCygOlo1.pri.ttmpnet GCF_000002315.6_GRCg6a.vs.GCA_009769625.1_bCygOlo1.pri.tnet

../software/kentutils/netSyntenic GCA_009769625.1_bCygOlo1.pri.qtmpnet GCA_009769625.1_bCygOlo1.pri.qnet


../software/kentut ils/netToAxt GCF_000002315.6_GRCg6a.vs.GCA_009769625.1_bCygOlo1.pri.tnet GCA_009769625.1_bCygOlo1.pri.prenet.chains GCF_000002315.6_GRCg6a_genomic.2bit GCA_009769625.1_bCygOlo1.pri_genomic.2bit GCA_009769625.1_bCygOlo1.pri.axt


##get the size information of sequences
for i in *_assembly_report.txt; do seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f7 | grep -v na`; if [[ ${#seqids} -eq 0 ]]; then seqids=`awk '$2=="assembled-molecule"' $i | grep -v non-nuclear | cut -f5 | grep -v na`; fi; (for j in $seqids; do echo -e "$j\t`grep -v ^# $i | grep $j | cut -f 9`"; done;) >`dirname $i`/`basename $i _assembly_report.txt`.sizes; done


grep -v ^# GCF_000002315.6_GRCg6a.vs.GCA_002798355.1_Ogye1.0.tnet | perl -lne 'if ($_=~/net (\S+) (\d+)/) { $tc=$1;$tcl=$2 } elsif ($_=~/^ fill (\d+) (\d+) (\S+) (\S) (\d+) (\d+)/) { $ts=$1; $tl=$2; $qc=$3; $strand = $4; $qs=$5; $ql=$6; $te=$ts+$tl; $qe=$qs+$ql; print "$tc\t$ts\t$te\t$qc\t$qs\t$qe\t$strand" }' >GCF_000002315.6_GRCg6a.vs.GCA_002798355.1_Ogye1.0.chainpairs.tab


Ensembl definition of Synteny
Synteny is the conserved order of aligned genomic blocks between species. It is calculated from the pairwise genome alignments created by Ensembl, when both species have a chromosome-level assembly.
The search is run in two phases:
1. We search for alignment blocks that are in the same order in the two genomes. Syntenic alignments that are closer than 200 kb are grouped into a synteny block.
2. Groups that are in synteny are linked, provided that no more than two non-syntenic groups are found between them and they are less than 3 Mb apart.
```
