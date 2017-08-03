#! /bin/bash

date
printf '\ndbs_automation: CHANGING SETTIGNGS...\n'

dbs_change_settings $1 type=WFA picardPath=$PICARD bowtie2Reference=/Users/pontushojer/data_analysis/Reference/Bowtie2Index/genome maxHandleMissMatches=4 barcodeMissmatch=2 mapqCutOff=20 port=5000 nexteraLayout=True

date
printf '\ndbs_automation: SETTINGS CHANGED\ndbs_automation: IDENTIFYING HANDLES...\n'
date > $1/logfiles/dbs_automation.txt
printf '\ndbs_automation: SETTINGS CHANGED\ndbs_automation: IDENTIFYING HANDLES...\n' >> $1/logfiles/dbs_automation.txt

dbs_handle_indetifier $2 $3 $1

date
printf '\ndbs_automation: HANDLES IDENTIFIED\ndbs_automation: CLUSTERING BARCODES...\n\n'
date >> $1/logfiles/dbs_automation.txt
printf '\ndbs_automation: HANDLES IDENTIFIED\ndbs_automation: CLUSTERING BARCODES...\n\n' >> $1/logfiles/dbs_automation.txt

dbs_barcode_clustrer $1

date
printf '\ndbs_automation: BARCODES CLUSTERED\ndbs_automation: MAPPING INSERTS...\n\n'
date >> $1/logfiles/dbs_automation.txt
printf '\ndbs_automation: BARCODES CLUSTERED\ndbs_automation: MAPPING INSERTS...\n\n' >> $1/logfiles/dbs_automation.txt

dbs_insert_mapper $1

date
printf '\ndbs_automation: INSERTS MAPPED\ndbs_automation: ANALYZING CLUSTERS...\n\n'
date >> $1/logfiles/dbs_automation.txt
printf '\ndbs_automation: INSERTS MAPPED\ndbs_automation: ANALYZING CLUSTERS...\n\n' >> $1/logfiles/dbs_automation.txt

dbs_cluster_analyzer $1

date
printf '\ndbs_automation: FINISHED.'
date >> $1/logfiles/dbs_automation.txt
printf '\ndbs_automation: FINISHED.' >> $1/logfiles/dbs_automation.txt
