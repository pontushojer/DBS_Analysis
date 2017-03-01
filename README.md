# DBS_Analysis

Software for analysis of Droplet Barcode Sequencing data.

## To install run:
```
    git clone https://github.com/elhb/DBS_Analysis.git
    cd DBS_Analysis
    git checkout structure
    pip install -r requirements.txt
    python setup.py build
    python setup.py install
```
**Note that:** *Some of the scripts make use of other non-python software you therefore need working installations of* **bowtie2**, **cd-hit-454** *and* **picard tools** *(version 1.114) to be able to run all commands.*
*It's easy to install all of this using conda, take a look at the analysis_automation.sh script.*
*Otherwise you can find the software here:*
- *bowtie2 can be found [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.8/)*
- *picard tools can be found [here](https://sourceforge.net/projects/picard/files/picard-tools/1.114/)*
- *cd-hit-454 can be found [here](https://github.com/weizhongli/cdhit/releases/download/V4.6.1/cd-hit-v4.6.1-2012-08-27.tgz)* 
**Note that:** *the OpenMP threading of cd-hit does not work for the conda installation.*

## To run the analysis of the 8 coriell individuals data:
This script is supplied to allow others to reproduce the analysis exactly as it was done for the publication.
```
    bash analysis_automation.sh
```
**Note that:** *This part is still under* **heavy** *development.*

## Standard Usage:

For each sample you want to analyze the following steps need to be performed ` <path to analysis> ` is the path to a folder where you want your analysis results and intermediary files to be stored.

#### 1. Initiate the analysis and define some values for adequate settings
*To review your settings or see what options to work with run the following command:*
```
    dbs_change_settings <path to analysis> listsettings
```
*To set a value for a variable run the following command:*
```
    dbs_change_settings <path to analysis> <variable name>=<value>
```
**Example:**
*Remember that the code supplied below is only an example and the paths most probably need to be changed, eg. the file "reference_data/hla_a.fasta" is probably not around on your computer at this moment.*
```
    dbs_change_settings ~/my_analysis/sample_one
        mapqCutOff=20 \
        maxHandleMissMatches=2 \
        maxIndividualIndexMissMatches=1 \
        barcodeMissmatch=2 \
        port=5000 \
        targetRegionBed=software/DBS_Analysis/goodies/target_regions.bed \
        IndexReferenceTsv=software/DBS_Analysis/goodies/indexRef.tsv \
        bowtie2Reference=reference_data/hla_a.fasta \
        picardPath=software/picard-tools-1.114/ \
        known_hla_types=reference_data/ipd.imgt.hla_A_gen.fasta \
        minReadDepthForVariantCalling=5 \
        minPairsPerCluster=20;
```

#### 2. Supply the input fastq files and look for handle and droplet barcode sequences
```
    dbs_handle_indetifier <path to read one fastq> <path to read two fastq> <path to analysis>;
```
**Note that:** *If you use the keyword *`SKIP`* instead of *`<path to read one fastq>`* files will not be added to the database but the handle and dbs identification will be run for all files already in the database:*
```
    dbs_handle_indetifier SKIP SKIP <path to analysis>;
```
**Note that:** *If you use the keyword *`justadd`* at the end of the line files only be added to the database and no handle and dbs identification will be performed:*
```
    dbs_handle_indetifier <path to read one fastq> <path to read two fastq> <path to analysis> justadd;
```
#### 3. Cluster the barcode sequences and group the reads based on this information
```
    dbs_barcode_clustrer <path to analysis>
```

#### 4. Map the reads to reference sequence while keeping the connections to barcode clusters in a database
```
    dbs_insert_mapper <path to analysis>
```
#### 5. Analyze each barcode cluster in parallel
```
    dbs_cluster_analyzer <path to analysis>
```
**Note that:** *It is sometimes usefull to be able to rerun this step from scratch by supplying the argument reanalyze=True:*
```
    dbs_cluster_analyzer <path to analysis> reanalyze=True
```
**Note that:** *This step sometimes freezez before finishing (eg. at 60% progress) this usually indicates that something is funky with the pysam installation, if this is the case try to build pysam from source, potentially also try switching virtualenv to conda this worked for me at several occasions.*

#### 6. Identify the alleles present within the data
```
    dbs_find_alleles <path to analysis>
```
**Note that:** *You can also supply the optional keyword graph to skip some of the computiationally intenst steps and just make some plots:*
```
    dbs_find_alleles <path to analysis> [graph]
```

#### 7. To view the results in a web browser
```
    dbs_hla_server <path to analysis>
```
## Other Usefull Commands:
#### 1. In case the IGV.js web browser vizualisation is slow it might help to dowsample the bamfiles
```
    dbs_downsample_bam <path to analysis>
```

#### 2. To add a barcode cluster specific tag to each read
this adds a barcode cluster specififc tag ("bc") to each read in the `<path to analysis>/data/mappedInserts.bam` creating a new file named `<path to analysis>/data/mappedInserts.tagged.bam`
```
    dbs_tag_bam <path to analysis>
```

#### 3. To view logfiles created during runtime
```
    viewLogfiles <path to analysis>/logfiles
```