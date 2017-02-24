# DBS_Analysis

Software for analysis of Droplet Barcode Sequencing data.

## To install run:
```
    git clone https://github.com/elhb/DBS_Analysis.git
    cd DBS_Analysis
    python setup.py build
    python setup.py install
```

## Standard Usage:

For each sample you want to analyze the following steps need to be performed ` <path to analysis> ` is the path to a folder where you want your analysis results and intermediary files to be stored.

##### 1. Initiate the analysis and define some values for adequate settings
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

##### 2. Supply the input fastq files and look for handle and droplet barcode sequences
```
    dbs_handle_indetifier <path to read one fastq> <path to read two fastq> <path to analysis>;
```

##### 3. Cluster the barcode sequences and group the reads based on this information
```
    dbs_barcode_clustrer <path to analysis>
```

##### 4. Map the reads to reference sequence while keeping the connections to barcode clusters in a database
```
    dbs_insert_mapper <path to analysis>
```
##### 5. Analyze each barcode cluster in parallel
```
    dbs_cluster_analyzer <path to analysis>
```
**Note that:** *It is sometimes usefull to be able to rerun this step from scratch by supplying the argument reanalyze=True:*
```
    dbs_cluster_analyzer <path to analysis> reanalyze=True
```

##### 6. Identify the alleles present within the data
```
    dbs_find_alleles <path to analysis>
```
**Note that:** *You can also supply the optional keyword graph to skip some of the computiationally intenst steps and just make some plots:*
```
    dbs_find_alleles <path to analysis> [graph]
```

##### 7. To view the results in a web browser
```
    dbs_hla_server <path to analysis>
```
## Other Usefull Commands:
##### 1. In case the IGV.js web browser vizualisation is slow it might help to dowsample the bamfiles
```
    dbs_downsample_bam <path to analysis>
```

##### 2. To add a barcode cluster specific tag to each read
this adds a barcode cluster specififc tag ("bc") to each read in the `<path to analysis>/data/mappedInserts.bam` creating a new file named `<path to analysis>/data/mappedInserts.tagged.bam`
```
    dbs_tag_bam <path to analysis>
```

##### 3. To view logfiles created during runtime
```
    viewLogfiles <path to analysis>/logfiles
```

## To run the analysis of the 8 coriell individuals data:
```
    still under construction
```
