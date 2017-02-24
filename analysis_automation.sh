#
# Prior to running this script you need a working install with executables in your path for the following software
#     virtualenv
#     bowtie2
#     cd-hit-454
#     git
#

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
# os specific block START ###########################################################################################################
#                                                                                                                                   #
# ubuntu version:                                                                                                                   #
# wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz"                                      #
# tar -xvzf sratoolkit.2.8.1-3-ubuntu64.tar.gz                                                                                      #
# sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-ubuntu64                                                                                        #
# wget "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip"                                  #
#                                                                                                                                   #
# centos version:                                                                                                                   #
# wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-centos_linux64.tar.gz"                                #
# tar -xvzf sratoolkit.2.8.1-3-centos_linux64.tar.gz                                                                                #
# sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-centos_linux64                                                                                  #
# wget "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip"                                  #
#                                                                                                                                   #
# OSX version:                                                                                                                      #
curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-mac64.tar.gz -o sratoolkit.2.8.1-3-mac64.tar.gz          #
tar -xvzf sratoolkit.2.8.1-3-mac64.tar.gz                                                                                           #
sratoolkitpath=$(pwd)/sratoolkit.2.8.1-3-mac64                                                                                             #
curl "https://netix.dl.sourceforge.net/project/picard/picard-tools/1.114/picard-tools-1.114.zip" -o picard-tools-1.114.zip          #
#                                                                                                                                   #
# os specific block END #############################################################################################################
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
for accession in SRR5277650 SRR5277651 SRR5277652 SRR5277653 SRR5277654 SRR5277655 SRR5277656 SRR5277657 SRR5277658;
    do $sratoolkitpath/bin/fastq-dump --split-files --gzip --outdir rawdata $accession; done
    
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
    echo -e "### to look at the results you can now run the following command:\n###\tdbs_hla_server " $(pwd)"/analysis_results/"$accession "\n### and then go to http://0.0.0.0:5000/ in your favourite web browser"
done;

deactivate