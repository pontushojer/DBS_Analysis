function main {
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m looking for a conda installation ... \033[0m";
    if [[ $(conda --version 2>&1) =~ conda ]];
    then
    
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m conda is installed creating an enviroment and installing dependencies\033[0m";
        setup_analysis_enviroment
        get_reference_data
        get_rawdata
        run_analysis
    
    else
    
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m conda not installed ...\033[0m"  
    
        install_conda
    
    fi
}



function install_conda {
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Trying to install miniconda.\033[0m"  
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m You will have to restart the analysis after this step (and also close you terminal window and open a new one)\033[0m"  
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m you are running linux, starting the ubuntu installer\033[0m"
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m you are running darwin, starting the OSX installer\033[0m"
        curl "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" -o Miniconda3-latest-MacOSX-x86_64.sh
        bash Miniconda3-latest-MacOSX-x86_64.sh
    else
        echo "Your operating system is not supported try running an ubuntu or OSX machine."
    fi
}


function setup_analysis_enviroment {
    
    STARTPATH=$(pwd)
    
    # create env
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m creating the virtual enviroment dbs_analysis\033[0m"
    conda create -n dbs_analysis python=2.7 --yes
    
    # start env
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m activating the virtual enviroment dbs_analysis (you can do this later by running: source activate dbs_analysis)\033[0m"
    source activate dbs_analysis
    
    # add channels and install dependencies
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m installing dependencies using conda\033[0m"
    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda
    conda install bowtie2  cd-hit cython flask numpy biopython Flask matplotlib pysam cutadapt wget --yes

    # create a directory to hold all the files needed for the analysis
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m create a directory ("$(pwd)"/analysis_automation) to hold all the files needed for the analysis\033[0m"
    mkdir -p analysis_automation
    cd analysis_automation
    
    # make a sub-directory that will hold the software we are going to use (apart from the dependencies installed earlier, this part should be more stable)
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Downloading and installing software\033[0m"
    mkdir software
    cd software
    
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Downloading and extracting the SRA toolkit\033[0m"
    
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
    
        # ubuntu version:
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Running ubuntu version\033[0m"
        
        # Download  and extract the sra toolkit
        curl --insecure "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz" -o sratoolkit.2.8.1-3-ubuntu64.tar.gz
        tar -xzf sratoolkit.2.8.1-3-ubuntu64.tar.gz
        sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-ubuntu64    
    elif [[ "$OSTYPE" == "darwin"* ]]; then
    
        # OSX version
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Running OSX version\033[0m"
    
        # Download  and extract the sra toolkit
        curl --insecure "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-mac64.tar.gz" -o sratoolkit.2.8.1-3-mac64.tar.gz
        tar -xvzf sratoolkit.2.8.1-3-mac64.tar.gz
        sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-mac64
    else
            echo "Your OS is not supported, aborting"; exit
    fi
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Downloading and extracting picard tools\033[0m"    
    # Download picardtools version 1.114 that was used in the publication
    curl --insecure -L "https://sourceforge.net/projects/picard/files/picard-tools/1.114/picard-tools-1.114.zip/download" -o picard-tools-1.114.zip
    # extract the content of the just downloaded picardtools
    unzip -q picard-tools-1.114.zip 
    
    # Now download and install our custom software that are specific for the droplet barcoding data
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Downloading and installing the DBS_Analysis software\033[0m"
    git clone https://github.com/elhb/DBS_Analysis.git
    cd DBS_Analysis/
    python setup.py build
    python setup.py install
    
    source deactivate dbs_analysis
    
    cd $STARTPATH
    
}


function get_reference_data {

    STARTPATH=$(pwd)
    mkdir -p analysis_automation
    cd analysis_automation

    # Getting reference data    
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Downloading and formating reference data\033[0m"
    mkdir reference_data
    curl --insecure "https://genome.ucsc.edu/cgi-bin/hgc?hgsid=580682449_q9IyW6Xnd8wjJrgFDMoAM9ORiftC&g=htcGetDna2&table=&i=mixed&o=29909823&l=29909823&r=29913852&getDnaPos=chr6%3A29907000-29917000&db=hg19&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA" | awk '/^>/{print $0} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/hla_a.fasta
    curl --insecure "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/A_gen.fasta" | awk '/^>/{print a;print $0;a=""} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/ipd.imgt.hla_A_gen.fasta
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Building a bowtie2 reference index\033[0m"
    source activate dbs_analysis
    bowtie2-build reference_data/hla_a.fasta reference_data/hla_a.fasta    
    source deactivate dbs_analysis
    cd $STARTPATH

}


function get_rawdata {

    STARTPATH=$(pwd)
    mkdir -p analysis_automation
    cd analysis_automation

    # Getting the raw sequence data
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Fetching the raw reads fastq files from SRA ... THIS MIGHT TAKE A WHILE!\033[0m"
    mkdir rawdata
    for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
        do $sratoolkitpath/bin/fastq-dump --split-files --gzip --outdir rawdata $accession; done    

    cd $STARTPATH
}


function run_analysis {
    # Now its time to run the analysis of our downloaded data,
    STARTPATH=$(pwd)
    mkdir -p analysis_automation
    cd analysis_automation
    source activate dbs_analysis
    # from here on the script also describes exactly how the analysis was made for the publication
    echo ""
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Running the analysis\033[0m"
    mkdir analysis_results
    
    # this might be needed on some headless systems to make the matplotlib plots work
    # export MPLBACKEND="agg" 
    
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m trimming the raw reads and create each analysis folder with appropriate settings\033[0m"
    # H2 CGATGAACTCGTGAAGCTAG should be trimmed from R2
    # Individual ID handle CTGAGTCCCGGTGGGTGCGTNNNNNNATTGAAGCCTGCCGTCGGAGACTAACGATGAACTCGTGAAGCTAG should also be trimmed from R2
    # H3 CGCTATCATAACGGCCATCT should be trimmed from R1
    for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
    do
        r1=$(echo $accession|awk '{print "rawdata/"$1"_1.fastq.gz"}');
        r2=$(echo $accession|awk '{print "rawdata/"$1"_2.fastq.gz"}');
        cutadapt -a CGCTATCATAACGGCCATCT -A CGATGAACTCGTGAAGCTAG -A CTGAGTCCCGGTGGGTGCGTNNNNNNATTGAAGCCTGCCGTCGGAGACTAACGATGAACTCGTGAAGCTAG -o $r1.trimmed.fq -p $r2.trimmed.fq $r1 $r2 >> rawdata/$accession.cutadapt_log.txt 2>&1 &
    
        dbs_change_settings analysis_results/$accession\
            mapqCutOff=20 \
            maxHandleMissMatches=2 \
            maxIndividualIndexMissMatches=1 \
            barcodeMissmatch=2 \
            port=5000 \
            targetRegionBed=$(pwd)/software/DBS_Analysis/goodies/target_regions.bed \
            IndexReferenceTsv=$(pwd)/software/DBS_Analysis/goodies/indexRef.tsv \
            bowtie2Reference=$(pwd)/reference_data/hla_a.fasta \
            picardPath=$(pwd)/software/picard-tools-1.114/ \
            known_hla_types=$(pwd)/reference_data/ipd.imgt.hla_A_gen.fasta \
            minReadDepthForVariantCalling=5 \
            minPairsPerCluster=20;
    done
    
    wait
    
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Compressing the trimmed data\033[0m"
    for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
    do
        r1=$(echo $accession|awk '{print "rawdata/"$1"_1.fastq.gz"}');
        r2=$(echo $accession|awk '{print "rawdata/"$1"_2.fastq.gz"}');
        gzip -v9  $r1.trimmed.fq &
        gzip -v9  $r2.trimmed.fq &
    done
    
    wait
    
    echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m For each sample run the analysis scripts\033[0m"
    for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
    do
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m Now processing sample $accession\033[0m"
        r1=$(echo $accession|awk '{print "rawdata/"$1"_1.fastq.gz"}');
        r2=$(echo $accession|awk '{print "rawdata/"$1"_2.fastq.gz"}');
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: searching for handles and barcode information\033[0m"
        dbs_handle_indetifier $r1.trimmed.fq.gz $r2.trimmed.fq.gz analysis_results/$accession;
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: cluster the barcode sequences and group the reads based on this information\033[0m"
        dbs_barcode_clustrer analysis_results/$accession;
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: mapping the reads to reference sequence\033[0m"
        dbs_insert_mapper analysis_results/$accession;
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: compressing intermediary files\033[0m"
        gzip -v9 analysis_results/$accession/data/inserts.singlets.fastq &
        gzip -v9 analysis_results/$accession/data/inserts.r1.fastq &
        gzip -v9 analysis_results/$accession/data/inserts.r2.fastq &
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: analyze each barcode cluster in parallel\033[0m"
        dbs_cluster_analyzer analysis_results/$accession;
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: identify the alleles present within the data\033[0m"
        dbs_find_alleles analysis_results/$accession >> analysis_results/$accession.dbs_find_alleles_output.txt 2>&1;
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession :: build some graphs\033[0m"
        dbs_find_alleles analysis_results/$accession graph > DELETETHISFILE.txt
        rm DELETETHISFILE.txt
        sleep 1
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m sample $accession is now finished!\033[0m"
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m to look at the results you can now run the following commands:\033[0m"
        echo -e "\033[1;93m##INFO##\033[0m ---\tsource activate dbs_analysis\033[0m"
        echo -e "\033[1;93m##INFO##\033[0m ---\tdbs_hla_server " $(pwd)"/analysis_results/"$accession "\033[0m"
        echo -e "\033[1;93m##INFO##\033[0m --- \033[0;34m and then go to \033[0m http://0.0.0.0:5000/ \033[0;34m in your favourite web browser\033[0m"
        #
        # Note that if the IGV.js vizualisation is slow or hangs it might be due to that there are simply to many reads in a small region so that the IGV.js downsampling can't cope with it
        #     if this is the case you could run the dbs_downsample bam command as below as a last step to only make the webinterface nicer
        #     note that after this step has been done you should not rerun any analysissteps without first manually restoring the orirginal bamfiles in the data folder
        # dbs_downsample_bam analysis_results/$accession
    done;
    source deactivate dbs_analysis

    cd $STARTPATH
}



main