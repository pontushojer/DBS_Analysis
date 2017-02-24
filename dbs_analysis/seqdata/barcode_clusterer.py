class BarcodeClusterer(object):
    
    def __init__(self, analysisfolder):
        
        self.analysisfolder = analysisfolder
        self.logfile = self.analysisfolder.logfile

    def generateBarcodeFastq(self, ):
        """ function that loads the barcode sequnces found from the database and creates a fasta file with these sequence formated and rady for running the clustering
        """

        #
        # imports
        #
        import operator
        from dbs_analysis import misc
        if self.analysisfolder.settings.type == 'HLA':  from dbs_analysis.sequences import HLA_DBS as DBS
        elif self.analysisfolder.settings.type == 'WFA':from dbs_analysis.sequences import WFA_DBS as DBS

        if self.logfile: self.logfile.write('Generating barcode fastq ...\n')

        #
        # setting initial values
        #
        uniqBarcodeSequences = {}
        temporaryDict = {}
        barcodeCounter = 0
        totalReadPairCounter = 0
        qualities = {}
        readPairHasBarcodeCounter = 0
        base_frequencies = [{'A':0,'T':0,'G':0,'C':0} for i in xrange(len(DBS))]

        if self.logfile: self.logfile.write('Loading read pairs ...\n')
        progress = misc.Progress(self.analysisfolder.results.totalReadCount, logfile=self.logfile, unit='reads-loaded-from-db', mem=True)
        with progress:
            self.analysisfolder.database.getConnection()
            for pairid,barcodeSequence,qual in self.analysisfolder.database.c.execute('SELECT id, dbsSeq, dbsQual FROM reads'):
                pairid = int(pairid)
                if barcodeSequence:
                    readPairHasBarcodeCounter += 1
                    qualities[pairid] = qual
                    try:            uniqBarcodeSequences[barcodeSequence].append(pairid)
                    except KeyError:uniqBarcodeSequences[barcodeSequence] = [pairid]
                    for i in xrange(len(barcodeSequence)): base_frequencies[i][barcodeSequence[i]] += 1
                progress.update()
        if self.logfile: self.logfile.write('Done.\n')
        self.analysisfolder.results.setResult('uniqueBarcodeSequences',len(uniqBarcodeSequences))
        if self.logfile: self.logfile.write(str(self.analysisfolder.results.uniqueBarcodeSequences)+' uniq barcode sequences found within the read pair population.\n')
        
        #
        # print base frequenzies for raw barcodes to file
        #
        with open(self.analysisfolder.dataPath+'/rawBarcodeBaseFreq.dict.txt','w') as outfile: outfile.write(str(base_frequencies))

        self.id2seq = {}
        if self.logfile: self.logfile.write('Sorting the barcodes by number of reads/sequence.\n')
        if self.logfile: self.logfile.write('Building the sorting dictionary ...\n')
        for barcode, idList in uniqBarcodeSequences.iteritems():
            try:		temporaryDict[len(idList)].append(barcode)
            except KeyError:temporaryDict[len(idList)] = [barcode]
            for idnumber in idList: self.id2seq[idnumber] = barcode

        if self.logfile: self.logfile.write('Creating output ... \n')
        barcodeFastqFile = open(self.analysisfolder.dbsfastq,'w')
        progress = misc.Progress(readPairHasBarcodeCounter, logfile=self.logfile, unit='reads-to-fastq', mem=True)
        with progress:
            for count, barcodes in sorted(temporaryDict.iteritems(), key=operator.itemgetter(0))[::-1]:
                for barcode in barcodes:
                    barcodeCounter += 1
                    readPairCounter = 0
                    for readPairId in uniqBarcodeSequences[barcode]:
                        readPairCounter += 1
                        totalReadPairCounter += 1
                        barcodeFastqFile.write('@'+str(readPairId)+' bc='+str(barcodeCounter)+' rp='+str(readPairCounter)+' bctrp='+str(count)+'\n'+barcode+'\n+\n'+qualities[readPairId]+'\n')
                        progress.update()
        barcodeFastqFile.close()	
        
        return readPairHasBarcodeCounter

    def runBarcodeClusteringPrograms(self, ):
        """ runs the actual clustering programs as subprocesses
        """

        #
        # imports
        #
        import subprocess
        import time

        clusteringProgramsLogfile = open(self.analysisfolder.logpath+'/'+time.strftime("%y%m%d-%H:%M:%S",time.localtime())+'_cdHitBarcodeClustering.log.txt','w')

        # cluster raw barcodes
        import platform
        import re
        if re.search('Ubuntu',platform.platform()):
            cdhit454name = 'cdhit-454'
        else:
            cdhit454name = 'cd-hit-454'
        command = [cdhit454name,
                  '-i',self.analysisfolder.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq',
                  '-o',self.analysisfolder.dataPath+'/clusteredBarcodeSequences',
                  '-g','1',      # mode
                  '-n','3',      # wrodsize
                  '-M','0',      # memory limit
                  '-T',str(self.analysisfolder.settings.parallelProcesses),      # threads
                  '-gap','-100', # disallow gaps
                  '-c',str(1.000-float(self.analysisfolder.settings.barcodeMissmatch)/float(self.analysisfolder.settings.barcodeLength)) # identity cutoff
                  ]
        if self.logfile: self.logfile.write('Starting command: '+' '.join(command)+'\n')
        cdhit = subprocess.Popen(command,stdout=clusteringProgramsLogfile,stderr=subprocess.PIPE )
        errdata = cdhit.communicate()
        if cdhit.returncode != 0:
            print 'cmd: '+' '.join( command )
            print 'cd-hit view Error code', cdhit.returncode, errdata
            sys.exit()

        # Build consensus sequences for read pair clusters
        # command = ['cdhit-cluster-consensus',
        #     self.analysisfolder.dataPath+'/clusteredBarcodeSequences.clstr',
        #     self.analysisfolder.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq',
        #     self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus',
        #     self.analysisfolder.dataPath+'/clusteredBarcodeSequences.aligned'
        #     ]
        # if self.logfile: self.logfile.write('Starting command: '+' '.join(command)+'\n')
        # ccc = subprocess.Popen(command,stdout=clusteringProgramsLogfile,stderr=subprocess.PIPE )
        # errdata = ccc.communicate()
        # if ccc.returncode != 0:
        #     print 'cmd: '+' '.join( command )
        #     print 'cdhit-cluster-consensus view Error code', ccc.returncode, errdata
        #     sys.exit()
        clusteringProgramsLogfile.close()

        return 0

    def parseBarcodeClusteringOutput(self, readPairsHasBarcode):
        """ parse the output from the clustering programs to find what reads ids have beeen clustered together
        """

        from dbs_analysis import misc
        
        #
        # inititate variables
        #
        totalClusterCount = 0
        barcodeClusters = {}
        singletonClusters = {}
        nonSingletonClusters = {}

        #
        # open file connections
        #
        # consensusFile = open(self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus.fastq')
        clstrFile = open(self.analysisfolder.dataPath+'/clusteredBarcodeSequences.clstr')

        #
        # load cluster ids and consensus sequences
        #
        if self.logfile: self.logfile.write('\nLoading barcode clusters and a barcode consesnsus sequences for each cluster ...\n')
        # while True:
        #     header = consensusFile.readline().rstrip()
        #     barcodeSequence = consensusFile.readline().rstrip()
        #     junk = consensusFile.readline().rstrip()
        #     barcodeQuality = consensusFile.readline().rstrip()
        #     if header == '': break
        #     totalClusterCount += 1
        #     header = header.split('_cluster_')
        #     clusterId = int(header[1].split(' ')[0])
        #     if header[0][:2] == '@s':
        #         singletonClusters[clusterId] = {'clusterReadCount':1,'readPairs':[],'identities':[],'clusterBarcodeSequence':barcodeSequence,'clusterBarcodeQuality':barcodeQuality}
        #         barcodeClusters[clusterId] = singletonClusters[clusterId]
        #     elif header[0][:2] == '@c':
        #         nonSingletonClusters[clusterId] = {'clusterReadCount':int(header[1].split(' ')[2]),'readPairs':[],'identities':[],'clusterBarcodeSequence':barcodeSequence,'clusterBarcodeQuality':barcodeQuality}
        #         barcodeClusters[clusterId] = nonSingletonClusters[clusterId]
        #     else: raise ValueError
        # self.analysisfolder.results.setResult('barcodeClusterCount',totalClusterCount)
        # self.analysisfolder.results.setResult('singeltonBarcodeClusters',len(singletonClusters))
        if self.logfile: self.logfile.write('A total of '+str(totalClusterCount)+' clusters of barcode sequences were loaded into memory.\n')

        #
        # Load what readpairs are in each cluster
        #
        if self.logfile: self.logfile.write('\nLoading read pair to barcode cluster connections ...\n')
        progress = misc.Progress(readPairsHasBarcode, logfile=self.logfile, unit='reads-loaded', mem=True)
        # tmp_totalClusterCount = 0
        barcodeCluster = None
        with progress:
            for line in clstrFile:
                line = line.rstrip()
                if line[0] == '>':
                    if barcodeCluster:
                        # print '#',clusterId,barcodeCluster['clusterReadCount'],barcodeCluster['clusterBarcodeSequence'],(clusterId in nonSingletonClusters),(clusterId in singletonClusters),len(barcodeCluster['readPairs']),len(barcodeCluster['identities'])
                        # assert barcodeCluster['clusterReadCount'] == barcodeClusters[clusterId]['clusterReadCount'],'readcounts dont match'
                        # if barcodeCluster['clusterBarcodeSequence'] != barcodeClusters[clusterId]['clusterBarcodeSequence']:print'NNNNNNNNNNNNNNNNNNNN';print barcodeCluster['clusterBarcodeSequence'];print barcodeClusters[clusterId]['clusterBarcodeSequence']
                        barcodeClusters[clusterId] = barcodeCluster
                        if barcodeCluster['clusterReadCount'] >1: nonSingletonClusters[clusterId] = barcodeCluster
                        else: singletonClusters[clusterId] = barcodeCluster
                    # tmp_totalClusterCount += 1
                    clusterId = int(line.split(' ')[1])
                    clusterReadCount = 0
                    barcodeCluster = {'clusterReadCount':0,'readPairs':[],'identities':[],'clusterBarcodeSequence':None,'clusterBarcodeQuality':None}
                    continue
                elif line[0] == '0':
                    readId = line.split('>')[1].split('.')[0]
                    identity = 'seed'
                    assert line.split(' ')[-1] == '*', 'Error in file format of clstr file'
                    barcodeCluster['clusterBarcodeSequence'] = self.id2seq[int(readId)]
                else:
                    readId = line.split('>')[1].split('.')[0]
                    identity = float(line.split(' ')[-1].split('/')[-1].split('%')[0])
                # barcodeClusters[clusterId]['readPairs'].append(readId)
                # barcodeClusters[clusterId]['identities'].append(identity)
                # clusterReadCount += 1
                barcodeCluster['readPairs'].append(int(readId))
                barcodeCluster['identities'].append(identity)
                barcodeCluster['clusterReadCount'] += 1
                progress.update()
        totalClusterCount = len(barcodeClusters)
        self.analysisfolder.results.setResult('barcodeClusterCount',len(barcodeClusters))
        self.analysisfolder.results.setResult('singeltonBarcodeClusters',len(singletonClusters))
        if self.logfile: self.logfile.write('All read pair to barcode cluster connections loaded.\n')

        return barcodeClusters

    def addBarcodeClusterInfoToDatabase(self, barcodeClusters):
        """ adds the cluster info to the database both one entry for each cluster but also updates each readpair entry
        """

        from dbs_analysis import misc
        from dbs_analysis import metadata

        #
        # set initial values
        #
        tmpUpdateValues = {}
        updateValues = []
        updateChunks = []
        updateChunkSize = 10000
        addValues = [(None,0,'[]','[]','NoneCluster',None,None,None)]

        #
        # drop old data and create table for new data
        #
        if self.logfile: self.logfile.write('Create barcodeClusters table (and drop old one if needed) ...\n')
        self.analysisfolder.database.getConnection()
        self.analysisfolder.database.c.execute("DROP TABLE IF EXISTS barcodeClusters")
        self.analysisfolder.database.c.execute('''CREATE TABLE barcodeClusters (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,PRIMARY KEY (clusterId))''')
        if self.logfile: self.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()

        #
        # Convert the dictionary to be able to add info to database
        #
        progress = misc.Progress(len(barcodeClusters), logfile=self.logfile, unit='clusters', mem=True)
        if self.logfile: self.logfile.write('Converting the data ... \n')
        for clusterId,data in barcodeClusters.iteritems():
            for readPairId in data['readPairs']: tmpUpdateValues[readPairId]=clusterId
            addValues.append( (clusterId,data['clusterReadCount'],str(data['readPairs']),str(data['identities']),data['clusterBarcodeSequence'],data['clusterBarcodeQuality'],None,None) )
            progress.update()
        for readPairId in sorted(tmpUpdateValues.keys()):
            updateValues.append( (int(tmpUpdateValues[readPairId]),int(readPairId)) )
            if len(updateValues) == updateChunkSize:
                updateChunks.append(updateValues)
                updateValues = []
        updateChunks.append(updateValues)

        #
        # Add the barcodeClusters
        #	
        if self.logfile: self.logfile.write('Adding cluster info to database ... \n')
        self.analysisfolder.database.getConnection()	
        self.analysisfolder.database.c.executemany('INSERT INTO barcodeClusters VALUES (?,?,?,?,?,?,?,?)', addValues)
        if self.logfile: self.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()

        #
        # Update the reads table
        #
        if self.logfile: self.logfile.write('Updating read pair info in the database ... \n')
        progress = misc.Progress(len(tmpUpdateValues), logfile=self.logfile, unit='reads-updated', mem=True)
        with progress:
            for updateValues in updateChunks:
                self.analysisfolder.database.getConnection()	
                self.analysisfolder.database.c.executemany('UPDATE reads SET clusterId=? WHERE id=?', updateValues)
                self.analysisfolder.database.commitAndClose()
                for i in xrange(len(updateValues)): progress.update()
            
        return 0

    def clusterBarcodeSequences(self):
        """ runns all the other functions in the right order
        """

        #
        # imports
        #
        import subprocess

        #
        # generate a fastq file with barcode sequences
        #
        readPairsHasBarcode = self.generateBarcodeFastq()
        self.analysisfolder.results.setResult('readPairsHasBarcode',readPairsHasBarcode)

        #
        # Run clustering programs
        #
        self.runBarcodeClusteringPrograms()

        #
        # Parse the output
        #
        barcodeClusters = self.parseBarcodeClusteringOutput(readPairsHasBarcode)

        #
        # compress files
        #
        a = [
              self.analysisfolder.dataPath+'/clusteredBarcodeSequences',
              self.analysisfolder.dataPath+'/clusteredBarcodeSequences.clstr',
              self.analysisfolder.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq',
              #self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus.fastq',
              #self.analysisfolder.dataPath+'/clusteredBarcodeSequences.aligned',
              ]
        processes = [] 
        if self.logfile: temp = self.logfile
        else: temp = sys.stdout
        for fileName in a:
            if self.logfile: self.logfile.write('Starting process:'+' gzip -v9 '+fileName+'\n')
            p=subprocess.Popen(['gzip','-f','-v9',fileName],stdout=temp)
            processes.append(p)

        #
        # save clustering info to database
        #
        self.addBarcodeClusterInfoToDatabase(barcodeClusters)

        #
        # Save results
        #
        self.analysisfolder.results.saveToDb()

        #
        # wait for file compression to finish
        #
        if self.logfile: self.logfile.write('Waiting for the compression of clustering files to finish.\n')
        for p in processes: p.wait()

        return 0

    def getBarcodeClusterIds(self, shuffle=True,byMixedClusterReadCount=True):
        """ function for fetching barcode cluster info from the database
        NOTE: this function should maybe not be under this object but rather be moved to the database object in the future(?)
        """

        import random
        from dbs_analysis.seqdata import BarcodeCluster

        self.analysisfolder.database.getConnection()
        tmp = list(self.analysisfolder.database.c.execute('SELECT clusterId,clusterTotalReadCount FROM barcodeClusters').fetchall())
        clusterIds = []
        readCountbyID = {}
        for clusterId,clusterTotalReadCount in tmp:
            clusterIds.append(clusterId)
            readCountbyID[clusterId] = clusterTotalReadCount
        if shuffle: random.shuffle(clusterIds)
        self.analysisfolder.database.commitAndClose()

        if byMixedClusterReadCount:
            
            id2ReadCountDist = {}
            yielded = {}
            normal = []
            large = []
            huge = []
            humongous = []
            outofthescale = []
            
            self.analysisfolder.logfile.write('Sorting cluster ids to groups of varying read pair counts ...\n')
            for clusterId in clusterIds:
                if clusterId == None: continue
                clusterId = int(clusterId)
                cluster = BarcodeCluster(clusterId,self.analysisfolder)
                #cluster.loadClusterInfo()
                cluster.readPairCount = readCountbyID[cluster.id]
                id2ReadCountDist[cluster.id] = cluster.readPairCount
                yielded[cluster.id] = False
                if   cluster.readPairCount <  1000 and cluster.readPairCount >= 0: normal.append(cluster.id)
                elif cluster.readPairCount <  5000 and cluster.readPairCount >= 1000: large.append(cluster.id)
                elif cluster.readPairCount < 35000 and cluster.readPairCount >= 5000: huge.append(cluster.id)
                elif cluster.readPairCount <500000 and cluster.readPairCount >= 35000: humongous.append(cluster.id)
                else: outofthescale.append(cluster.id)
            self.analysisfolder.logfile.write('Done.\n')
            
            tmpCounter0 = 0
            yieldNow = None
            while False in yielded.values():
                if outofthescale and tmpCounter0 <= 0:
                    yieldNow = outofthescale[0]
                    outofthescale = outofthescale[1:]
                    yielded[yieldNow] = True
                    tmpCounter0 += 1
                    yield yieldNow
                elif humongous and tmpCounter0 <= 0:
                    yieldNow = humongous[0]
                    humongous = humongous[1:]
                    yielded[yieldNow] = True
                    tmpCounter0 += 1
                    yield yieldNow
                elif huge and tmpCounter0 <= 0:
                    yieldNow = huge[0]
                    huge = huge[1:]
                    yielded[yieldNow] = True
                    tmpCounter0 += 1
                    yield yieldNow
                elif large and tmpCounter0 <= 1:
                    yieldNow = large[0]
                    large = large[1:]
                    yielded[yieldNow] = True
                    tmpCounter0 += 1
                    yield yieldNow
                elif normal:
                    yieldNow = normal[0]
                    normal = normal[1:]
                    yielded[yieldNow] = True
                    tmpCounter0 += 1
                    yield yieldNow
                #SEAseqPipeLine.logfile.write( str(tmpCounter0)+'\t')
                #SEAseqPipeLine.logfile.write( str(yieldNow)+'\t')
                #SEAseqPipeLine.logfile.write( str(id2ReadCountDist[yieldNow])+'.\n')
                if tmpCounter0 == 6: tmpCounter0 = 0
                #if yielded.values().count(True) >= 100: break
        else:
            for clusterId in clusterIds:
                if clusterId[0] == None: continue
                try: yield int(clusterId[0])
                except TypeError: yield clusterId[0]