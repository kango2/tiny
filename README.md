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

**List of species considered for analyses**
|Species|Species Code|Common Name|NCBI Taxonomy ID|Citation|
|:---|:---|:---|---:|:---|
|Gallus gallus|CHICK|Chicken|9031|[Hillier, L.W. et al, 2004](https://doi.org/10.1038/nature03154)|
|Trachemys scripta elegans|TRASE|Red-eared slider turtle|31138|[Simison, W.B. et al, 2020](https://doi.org/10.1093/gbe/evaa063)|
|Corvus monedula|CORMO|Jackdaw|30423|[Weissensteiner, M.H. et al, 2020](https://doi.org/10.1038/s41467-020-17195-4)|
|Cygnus olor|CYGOL|Mute swan|8869|VGP & B10K: Embargoed|
|Naja naja |NAJNA|Indian cobra|35670|[Kushal Suryamohan, K. et al, 2020](https://doi.org/10.1038/s41588-019-0559-8)|
|Ornithorhynchus anatinus|ORNAN|Platypus|9258|[Zhou, Y., et al, 2021](https://doi.org/10.1038/s41586-020-03039-0)|
|Gopherus evgoodei|GOPEV|Goodes thornscrub tortoise|1825980|VGP: Embargoed|
|Dermochelys coriacea|DERCO|Leatherback sea turtle|27794|VGP: Embargoed|
|Chelonia mydas|CHEMY|Green sea turtle|8469|VGP: Embargoed.<br />[Also at DNAzoo](https://www.dnazoo.org/assemblies/Chelonia_mydas) from male blood. Draft genome: [Zhuo, W. et al, 2013](https://doi.org/10.1038/ng.2615)|
|Lacerta agilis|LACAG|Sand lizard|80427|VGP: Embargoed|
|Thamnophis elegans|THAEL|Western terrestrial garter snake|35005|VGP: Embargoed|
|Podarcis muralis|PODMU|Common wall lizard|64176|[Andrade, P. et al, 2019](https://doi.org/10.1073/pnas.1820320116)|
|Zootoca vivipara|ZOOVI|Common lizard|8524|[Yurchenko, A.A. et al, 2020](https://doi.org/10.1093/gbe/evaa161)|
|Crotalus viridis|CROVV|Prairie rattlesnake|8742|[Pasquesi, G.I.M. et al, 2018](https://doi.org/10.1038/s41467-018-05279-1)|
|Calypte anna|CALAN|Anna's hummingbird|9244|[Rhie, A. et al, 2021](https://doi.org/10.1101/2020.05.22.110833)|
|Falco rusticolus|FALRU|Gyrfalcon|120794|VGP & B10K: Embargoed|
|Aquila chrysaetos chrysaetos|AQUCH|Golden eagle|223781|Sanger 25G & VGP: Embargoed.<br />[Also at DNAzoo](https://www.dnazoo.org/assemblies/Aquila_chrysaetos) from female blood. Draft genome: [Bussche, R.A. et al, 2017](https://doi.org/10.3356/JRR-16-47.1)|
|Dromaius novaehollandiae|DRONO|Emu|8790|[DNAzoo](https://www.dnazoo.org/assemblies/Dromaius_novaehollandiae) from male blood. Draft genome: [Sackton, T.B. et al, 2019](https://doi.org/10.1126/science.aat7244)|
|Salvator merianae|SAMER|Argentine black and white tegu|96440|[DNAzoo](https://www.dnazoo.org/assemblies/Salvator_merianae) from male blood. Draft genome: [Roscito, J.G. et al, 2018](https://doi.org/10.1093/gigascience/giy141)|
|Python bivittatus|PYTBI|Burmese python|176946|[DNAzoo](https://www.dnazoo.org/assemblies/Python_bivittatus) from feamle blood. Draft genome: [Castoe, T.A. et al, 2013](https://doi.org/10.1073/pnas.1314475110)|
|Gymnogyps californianus|GYMCA|California condor|33616|[DNAzoo](https://www.dnazoo.org/assemblies/Gymnogyps_californianus) from blood. Draft genome: Unpublished|
|Alligator mississippiensis|ALLMI|American alligator|8496|[DNAzoo](https://www.dnazoo.org/assemblies/Alligator_mississippiensis) from male blood. Draft genome: [St John, J.A., et al, 2017](https://doi.org/10.1101/gr.213595.116)|
|Alligator sinensis|ALLSI|Chinese alligator|38654|[DNAzoo](https://www.dnazoo.org/assemblies/Alligator_sinensis) from female blood. Draft genome: [Wan, Q.H. et al, 2013](https://doi.org/10.1038/cr.2013.104)|
|Homo sapiens|HUMAN|Human|9606|[Lander, E. et al, 2001](https://doi.org/10.1038/35057062)|
|Sarcophilus harrisii|SARHA|Tasmanian devil|9305|Wellcome Sanger Institute|
|Branchiostoma floridae|BRAFL|Florida lancelet|7739|[Simakov, O. et al, 2020](https://doi.org/10.1038/s41559-020-1156-z)|
|Phascolarctos cinereus|PHACI|Koala|38626|[DNAzoo](https://www.dnazoo.org/assemblies/Phascolarctos_cinereus) from male heart. Draft genome: [Johnson, R.N. et al, 2018](https://doi.org/10.1038/s41588-018-0153-5)|
|Casuarius casuarius|CASCA|Southern cassowary|8787|[DNAzoo](https://www.dnazoo.org/assemblies/Casuarius_casuarius) from male blood. Draft genome: [Sackton, T.B. et al, 2019](https://doi.org/10.1126/science.aat7244)|
|Tympanuchus cupido|TYMCU|Greater prairie chicken|9004|[DNAzoo](https://www.dnazoo.org/assemblies/Tympanuchus_cupido) from male blood. Draft genome: Unpublished|
|Eopsaltria australis|EOPAU|Eastern yellow robin|44318|[DNAzoo](https://www.dnazoo.org/assemblies/Eopsaltria_australis) from female liver. Draft genome: [Gan, H.M. et al, 2019](https://doi.org/10.1093/gigascience/giz111)|
|Lichenostomus melanops cassidix|LIMCA|Helmeted honeyeater|1497555|[DNAzoo](https://www.dnazoo.org/assemblies/Lichenostomus_melanops_cassidix) Draft genome: Unpublished|
|Patagioenas fasciata|PATFA|Band-tailed pigeon|372321|[DNAzoo](https://www.dnazoo.org/assemblies/Patagioenas_fasciata) from 'frozen' male sample. Draft genome: [Murray, G.G.R. et al, 2017](https://doi.org/10.1126/science.aao0960)|
|Phalacrocorax auritus|PHAAI|Double-crested cormorant|56069|[DNAzoo](https://www.dnazoo.org/assemblies/Phalacrocorax_auritus) from blood. Draft genome: [Burga, A. et al, 2017](https://doi.org/10.1126/science.aal3345)|
|Rhea americana|RHEAM|Greater rhea|8797|[DNAzoo](https://www.dnazoo.org/assemblies/Rhea_americana) from female blood. Draft genome: [Sackton, T.B. et al, 2019](https://doi.org/10.1126/science.aat7244)|
|Strix occidentalis|STROC|Spotted owl|201991|[DNAzoo](https://www.dnazoo.org/assemblies/Strix_occidentalis) sample not described. Draft genome: [Hanna, Z.R. et al, 2017](https://doi.org/10.1093/gbe/evx158)|
|Struthio camelus|STRCA|Common ostrich|8801|[DNAzoo](https://www.dnazoo.org/assemblies/Struthio_camelus) from female blood. Draft genome: [Zhang, G. et al, 2014](https://doi.org/10.1126/science.1251385)|
|Intellagama lesueurii lesueurii|INTLE|Eastern water dragon|103694|[DNAzoo](https://www.dnazoo.org/assemblies/Intellagama_lesueurii_lesueurii) from liver. Draft genome: Unpublished|


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
