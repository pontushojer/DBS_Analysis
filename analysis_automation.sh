#
# Prior to running this script you need a working install with executables in your path for the following software
#     virtualenv
#     bowtie2
#     cd-hit-454
#     git
#
if [[ "$OSTYPE" == "linux-gnu" ]]; then
# on ubuntu this is really easy, just run:
#   sudo apt-get install virtualenv bowtie2 python-dev cd-hit git python-tk;
    echo "### Running ubuntu version, installing dependencies"
    sudo apt-get install virtualenv bowtie2 python-dev cd-hit git python-tk
elif [[ "$OSTYPE" == "darwin"* ]]; then
    now_at=$(pwd)
    echo "### Running OSX version, installing dependencies"
    # git is installed by default in Xcode and all will fail if xcode is not installed so hopefully you have Xcode (or at least the command-line developer tools)
    
    xcode-select --install
    if [[ $(brew --version) =~ Homebrew ]];
    then
        echo "### using Home brew to install some stuff.";
        brew install bowtie2
        brew install cd-hit
    else
        echo "### homebrew not found, will try to install by other method";

        if [[ $PATH =~ ~/bin ]];
        then
            echo "### ~/bin folder in PATH, nice";
        else
            echo "### ~/bin folder not in PATH";
            mkdir -p ~/bin
            echo -e "### Attemting to add \"export PATH=\$HOME/bin:\$PATH\" to ~/.bashrc if you don't use ~/.bashrc-file you might have to fix this manually";
            echo -e "export PATH=\$HOME/bin:\$PATH" >> ~/.bashrc
        fi
        
        # bowtie2
        if [[ $(bowtie2 --version) =~ "bowtie2-align-s version" ]];
        then
            echo "### bowtie2 already installed.";
        else
            echo "### bowtie2 not installed, installing";
            cd ~/bin
            curl "https://excellmedia.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.8/bowtie2-2.2.8-macos-x86_64.zip" -o bowtie2-2.2.8-macos-x86_64.zip
            unzip bowtie2-2.2.8-macos-x86_64.zip
            ln -s ~/bin/bowtie2-2.2.8/bowtie2 ~/bin2/
            ln -s ~/bin/bowtie2-2.2.8/bowtie2-build ~/bin/
        fi
    
        # cdhit
        if [[ $(cd-hit-454 --version) =~ "CD-HIT version" ]];
        then
            echo "### cd-hit-454 installed";
        else
            echo "### cd-hit-454 not installed, installing ... NOT WORKING DO MANUALLY";
            cd ~/bin
            git clone https://github.com/weizhongli/cdhit.git
            cd cdhit/
            make
            ln -s ~/bin2/cdhit/cd-hit-454 ~/bin/
            cd ..
        fi

    fi

    #virualenvn
    echo "### installing virtualenv";
    pip install virtualenv
    cd $now_at
elif [[ "$OSTYPE" == "freebsd"* ]]; then
    echo "OS freebsd is not supported, aborting"; exit
else
    echo "Your OS is not supported, aborting"; exit
fi
# on centos maybe like this? have not tested it ...
# sudo yum install python-devel python-setuptools python-pip; sudo pip install --upgrade pip; sudo pip install virtualenv

#
# Setting and setting up the analysis enviroment
#
echo "### Starting"
mkdir analysis_automation
cd analysis_automation

echo ""
echo "### Creating and starting virtual enviroment ..."
virtualenv --python=python2.7 analysis_automation_venv
source analysis_automation_venv/bin/activate
pip install --upgrade pip
echo "### virtual env with name analysis_automation_venv created and started."

echo ""
echo "### Downloading and installing software"
mkdir software
cd software

echo ""
echo "### Downloading SRA toolkit and picardtools"

if [[ "$OSTYPE" == "linux-gnu" ]]; then

    ### ubuntu version:
    echo "### Running ubuntu version"
    wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz"
    tar -xvzf sratoolkit.2.8.1-3-ubuntu64.tar.gz
    sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-ubuntu64
    wget "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip"

    ### centos version: MIGHT BE CRAP
    # echo "### Running centos version"
    # wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-centos_linux64.tar.gz"
    # tar -xvzf sratoolkit.2.8.1-3-centos_linux64.tar.gz
    # sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-centos_linux64
    # wget "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip"
    
elif [[ "$OSTYPE" == "darwin"* ]]; then

    # OSX version
    echo "### Running OSX version"
    curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-mac64.tar.gz -o sratoolkit.2.8.1-3-mac64.tar.gz
    tar -xvzf sratoolkit.2.8.1-3-mac64.tar.gz
    sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-mac64
    curl "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip" -o picard-tools-1.114.zip

elif [[ "$OSTYPE" == "cygwin" ]]; then
        # POSIX compatibility layer and Linux environment emulation for Windows
        echo "OS cygwin is not supported, aborting"; exit
elif [[ "$OSTYPE" == "msys" ]]; then
        # Lightweight shell and GNU utilities compiled for Windows (part of MinGW)
        echo "OS msys is not supported, aborting"; exit
elif [[ "$OSTYPE" == "win32" ]]; then
        # I'm not sure this can happen.
        echo "OS win32 is not supported, aborting"; exit
elif [[ "$OSTYPE" == "freebsd"* ]]; then
        pip install srapy
        echo "OS freebsd is not supported, aborting"; exit
else
        echo "Your OS is not supported, aborting"; exit
fi

unzip picard-tools-1.114.zip 

echo ""
echo "### Downloading and installing the DBS_Analysis software"
git clone https://github.com/elhb/DBS_Analysis.git
cd DBS_Analysis/
git checkout structure
pip install -r requirements.txt
python setup.py build
python setup.py install
cd ../..

#
# Getting reference data
#
echo ""
echo "### Downloading and formating reference data"
mkdir reference_data
# wget "https://genome.ucsc.edu/cgi-bin/hgc?hgsid=580682449_q9IyW6Xnd8wjJrgFDMoAM9ORiftC&g=htcGetDna2&table=&i=mixed&o=29909823&l=29909823&r=29913852&getDnaPos=chr6%3A29907000-29917000&db=hg19&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA" | awk '/^>/{print $0} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/hla_a.fasta 
# wget "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/A_gen.fasta" | awk '/^>/{print a;print $0;a=""} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/ipd.imgt.hla_A_gen.fasta
curl "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/A_gen.fasta" | awk '/^>/{print a;print $0;a=""} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/ipd.imgt.hla_A_gen.fasta
curl "https://genome.ucsc.edu/cgi-bin/hgc?hgsid=580682449_q9IyW6Xnd8wjJrgFDMoAM9ORiftC&g=htcGetDna2&table=&i=mixed&o=29909823&l=29909823&r=29913852&getDnaPos=chr6%3A29907000-29917000&db=hg19&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA" | awk '/^>/{print $0} /^[AGTCN]+$/{a=a $0;} END { print a }' > reference_data/hla_a.fasta 
echo "### Building a bowtie2 reference index"
bowtie2-build reference_data/hla_a.fasta reference_data/hla_a.fasta

#
# Getting the raw sequence data
#
echo ""
echo "### Fetching the raw reads fastq files from SRA ... THIS MIGHT TAKE A WHILE!"
mkdir rawdata

if [[ "$OSTYPE" == "freebsd"* ]]; then
        cd rawdata
        get-project-sras.py -p 376266
        for file in $(ls -lh *.sra | awk '{print $9}'); do mv -v $file* $(echo $file | awk '{print substr($0,0,10)".sra"}'); done
        # problem ... cannot convert .sra to fastq on freebsd :(
        cd ..
else
    for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
        do $sratoolkitpath/bin/fastq-dump --split-files --gzip --outdir rawdata $accession; done
fi
    
#
# run the analysis
#
echo ""
echo "### Running the analysis"
mkdir analysis_results

# this might be needed on some headless systems to make the matplotlib plots work
# export MPLBACKEND="agg" 

echo "### trimming the raw reads and create each analysis folder with appropriate settings"
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

echo "### Compressing the trimmed data"
for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
do
    r1=$(echo $accession|awk '{print "rawdata/"$1"_1.fastq.gz"}');
    r2=$(echo $accession|awk '{print "rawdata/"$1"_2.fastq.gz"}');
    gzip -v9  $r1.trimmed.fq &
    gzip -v9  $r2.trimmed.fq &
done

wait

echo "### For each sample run the analysis scripts"
for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
do
    echo "### Now processing sample $accession"
    r1=$(echo $accession|awk '{print "rawdata/"$1"_1.fastq.gz"}');
    r2=$(echo $accession|awk '{print "rawdata/"$1"_2.fastq.gz"}');
    echo "### sample $accession :: searching for handles and barcode information"
    dbs_handle_indetifier $r1.trimmed.fq.gz $r2.trimmed.fq.gz analysis_results/$accession;
    sleep 1
    echo "### sample $accession :: cluster the barcode sequences and group the reads based on this information"
    dbs_barcode_clustrer analysis_results/$accession;
    sleep 1
    echo "### sample $accession :: mapping the reads to reference sequence"
    dbs_insert_mapper analysis_results/$accession;
    sleep 1
    echo "### sample $accession :: compressing intermediary files"
    gzip -v9 analysis_results/$accession/data/inserts.singlets.fastq &
    gzip -v9 analysis_results/$accession/data/inserts.r1.fastq &
    gzip -v9 analysis_results/$accession/data/inserts.r2.fastq &
    sleep 1
    echo "### sample $accession :: analyze each barcode cluster in parallel"
    dbs_cluster_analyzer analysis_results/$accession;
    sleep 1
    echo "### sample $accession :: identify the alleles present within the data"
    dbs_find_alleles analysis_results/$accession >> analysis_results/$accession.dbs_find_alleles_output.txt 2>&1;
    sleep 1
    echo "### sample $accession :: build some graphs"
    dbs_find_alleles analysis_results/$accession graph
    sleep 1
    echo "### sample $accession is now finished!"
    echo -e "### to look at the results you can now run the following commands:"
    echo -e "###\tsource analysis_automation/analysis_automation_venv/bin/activate"
    echo -e "###\tdbs_hla_server " $(pwd)"/analysis_results/"$accession ""
    echo -e "### and then go to http://0.0.0.0:5000/ in your favourite web browser"
done;

deactivate