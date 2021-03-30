# Micro-chromosome conservation in vertebrates

Contributors: Hardip Patel, Paul Waters, Arthur Georges, Jenny Graves

## Data preparation

The process was automated using the `preparegenomes.sh` script. Basic steps were as follows.

1. Genome assemblies were downloaded from relevant sources (NCBI or DNAZoo). Check `metadata/species.txt` for details of species and download paths.
2. `.2bit`, `.capsule`, `.sizes` files were created.
3. Genome files were split into smaller regions (1Mb sequence size and ~5Mb of total sequence per file) without overlap. **This step allows for embarassingly parallel lastz alignments to capitalise on large HPC facilities.**

## One-way all-vs-all whole genome alignments using [LastZ](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#intro)

We used the following parameters for LastZ alignments:  
```K=2400 L=3000 Y=9400 H=2000 --ambiguous=iupac``` 

**List of species included in analyses**
|Species|Species Code|Common Name|NCBI Taxonomy ID|Citation|
|:---|:---|:---|---:|:---|
|Gallus gallus|CHICK|Chicken|9031||
|Trachemys scripta elegans|TRASE|Red-eared slider turtle|31138|https://doi.org/10.1093/gbe/evaa063|
|Gallus gallus|CHICKYO|Chicken Yeonsan Ogye|9031|
|Corvus monedula|CORMO|Jackdaw|30423|https://doi.org/10.1038/s41467-020-17195-4|
|Cygnus olor|CYGOL|Mute swan|8869|Genome10K??|
|Naja naja |NAJNA|Indian cobra|35670|https://doi.org/10.1038/s41588-019-0559-8|
|Ornithorhynchus anatinus|ORNAN|Platypus|9258||
|Gopherus evgoodei|GOPEV|Goodes thornscrub tortoise|1825980|Genome10K??|
|Dermochelys coriacea|DERCO|Leatherback sea turtle|27794||
|Chelonia mydas|CHEMY|Green sea turtle|8469||
|Lacerta agilis|LACAG|Sand lizard|80427||
|Thamnophis elegans|THAEL|Western terrestrial garter snake|35005||
|Podarcis muralis|PODMU|Common wall lizard|64176||
|Zootoca vivipara|ZOOVI|Common lizard|8524||
|Crotalus viridis|CROVV|Prairie rattlesnake|8742||
|Calypte anna|CALAN|Anna's hummingbird|9244||
|Falco rusticolus|FALRU|Gyrfalcon|120794||
|Aquila chrysaetos chrysaetos|AQUCH|Golden eagle|223781||
|Dromaius novaehollandiae|DRONO|Emu|8790||
|Salvator merianae|SAMER|Argentine black and white tegu|96440||
|Python bivittatus|PYTBI|Burmese python|176946||
|Gymnogyps californianus|GYMCA|California condor|33616||
|Alligator mississippiensis|ALLMI|American alligator|8496||
|Alligator sinensis|ALLSI|Chinese alligator|38654||
|Homo sapiens|HUMAN|Human|9606||
|Sarcophilus harrisii|SARHA|Tasmanian devil|9305||
|Branchiostoma floridae|BRAFL|Florida lancelet|7739||

**Typical workflow**

This workflow is automated using `getMAF.sh` script.

1. Fix LastZ alignments to assign genomic coordinates to sequence alignments as alignments were performed for 1Mb sub-sequences (`fixlastz.pl`).
2. Performing chaining of LastZ alignments (`axtChain`).
3. Sort chains (`chainSort`).
4. Generate prenet files from chain output (`chainPreNet`).
5. Perform netting of chains (`chainNet`).
6. Report nets of alignments against target genome (`netSyntenic`).
7. Report nets of alignments against query genome (`netSyntenic`).

Further steps that can be peformed but not performed as yet for this project. Commands are included in the `getMAF.sh` script, however, we have commented them out for now.

1. Convert nets to `.axt` format (`netToaxt`).
2. Sort `.axt` alignment files (`axtSort`).
3. Generate `.maf` format output for alignments (`axtToMaf`).


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
