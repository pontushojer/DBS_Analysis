class Database(object):
    
    def __init__(self, dbPath, analysisfolder):
        self.path = dbPath
        self.analysisfolder = analysisfolder

        # creates a lock for acces probably not really needed though it might be a nice feature in the future
        import multiprocessing
        import ctypes
        manager = multiprocessing.Manager()
        self.lock = manager.RLock()
        self.writeInProgress = manager.Value(ctypes.c_bool,False)

    def getConnection(self,):
        #
        # Import useful stuff
        #
        import sqlite3
        import sys

        #
        # Create database and set
        #
        #sys.stderr.write('connect 1\n')
        try:
            #sys.stderr.write('connect 2\n')
            self.conn = sqlite3.connect(self.path)
            #sys.stderr.write('connect 3\n')

        except sqlite3.OperationalError:
            #sys.stderr.write('connect ALTERNATICE\n')
            print 'ERROR: Trouble with the database, plase check your commandline.'
            sys.exit()
        #sys.stderr.write('connect 4\n')
        self.c = self.conn.cursor()
        #sys.stderr.write('connect 5\n')
    
    def commitAndClose(self,):
        #
        # commit changes and close connection
        #
        self.conn.commit()
        self.conn.close()
    
    def create(self,):
        """ creates the database holding all information used in the analysis """

        self.getConnection()

        #
        # Create tables
        #
        self.c.execute('''CREATE TABLE runs (startTime,command,commandLine,finishedSuccessfully,masterPid)''')
        self.c.execute('''CREATE TABLE fastqs (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength,PRIMARY KEY (filePairId))''');
        self.c.execute('''CREATE TABLE settings (variableName,defaultValue,value,setTime,PRIMARY KEY (variableName))''')
        self.c.execute('''CREATE TABLE results (resultName,defaultValue,value,setTime,PRIMARY KEY (resultName))''')
  
        self.commitAndClose()

        import os
        os.chmod(self.path, 0664)

    def addToRunsTable(self, startTime, command, commandLine, finishedSuccessfully, masterPid):
        
        self.getConnection()
        
        #
        # check if pid already in database
        #
        t = (masterPid,)
        data = self.c.execute('SELECT masterPid, startTime FROM runs WHERE masterPid=?',t).fetchall()        
        if data:
            for tmp1,tmp2 in data:

        #
        # if pid and startTime matches update the "finishedSuccessfully" entry
        #
                if tmp1 == masterPid and tmp2 == startTime:
                    values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
                    self.c.execute('UPDATE runs SET finishedSuccessfully=? WHERE masterPid=? AND startTime=?', (finishedSuccessfully,masterPid,startTime))
        
        #
        # if not in the database add a new row
        #
        else:
            values = (startTime, command, commandLine, finishedSuccessfully, masterPid)
            self.c.execute('INSERT INTO runs VALUES (?,?,?,?,?)', values)
        
        self.commitAndClose()
        
        return 0
    
    def addFastqs(self, fastq1, fastq2,logfile=None):
        
        #
        # Imports
        #
        import sys
        from dbs_analysis import misc
        import os
        
        #
        # open connection to database
        #
        self.getConnection()
        
        filePairId = None
        filePairIds = []
        
        #
        # check if any of the fastqs already in database
        #
        data = self.c.execute('SELECT filePairId,fastq1,fastq2 FROM fastqs').fetchall()
        if data:
            for filePair in data:
                filePairId = int(filePair[0])
                filePairIds.append(filePairId)
                for fastq in [fastq1, fastq2]:
                    if os.path.basename(fastq) in [os.path.basename(filePair[1]),os.path.basename(filePair[2])]:
                        message = 'ERROR: '+fastq+' already in the database.\nExiting after error.'
                        print message
                        if logfile: logfile.write(message+'\n')
                        sys.exit(1)
        #
        # if not in the database add a new row
        #
        if logfile: logfile.write('Linking files to analysis path ... \n')
        os.symlink(fastq1, self.analysisfolder.rawdataPath+'/'+os.path.basename(fastq1))
        os.symlink(fastq2, self.analysisfolder.rawdataPath+'/'+os.path.basename(fastq2))
        if logfile: logfile.write('Getting readcount for file'+os.path.basename(fastq1)+' ... \n')
        readCount = misc.bufcount(fastq1)/4 #one read is four lines
        if logfile: logfile.write('...done. The file has '+str(readCount)+' reads.\n')
        addedToReadsTable = False#SEAseqPipeLine.startTimeStr
        minReadLength = 'NA'

        if filePairIds: filePairId = max(filePairIds)+1
        else: filePairId = 0
        values = (filePairId,self.analysisfolder.rawdataPath+'/'+os.path.basename(fastq1),self.analysisfolder.rawdataPath+'/'+os.path.basename(fastq2),readCount,addedToReadsTable,minReadLength)
        self.c.execute('INSERT INTO fastqs VALUES (?,?,?,?,?,?)', values)
        
        self.commitAndClose()
        
        return 0
   
    def addReads(self, readsToAdd):

        #
        # Imports
        #
        import sys
        
        with self.lock:
            #
            # open connection to database
            #
            self.getConnection()
            
            #
            # add the data in readsToAdd to the reads table
            #
    
            self.c.executemany('INSERT INTO reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', readsToAdd)
            
            self.commitAndClose()
        
        return 0

    def getFastqs(self,):
        #
        # Imports
        #
        import sys
        
        #
        # open connection to database
        #
        self.getConnection()
                
        #
        # get att data in fastqs table
        #
        filePairs = self.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength FROM fastqs').fetchall()
        
        self.commitAndClose()
        
        #return [[readCount,fastq1,fastq2] if (not addedToReadsTable) else None for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength in filePairs]
        return [[filePairId,readCount,fastq1,fastq2] for filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength in filePairs]

    def getAllReadPairs(self,):
        #
        # Imports
        #
        import sys
        from dbs_analysis import seqdata
        
        #
        # open connection to database
        #
        self.getConnection()
                
        #
        # get att data in fastqs table
        #
        readPairs = self.c.execute('SELECT id, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id FROM reads')
        
        while True:

            rows = readPairs.fetchmany()#size=readPairs.arraysize)

            if not rows: break

            for row in rows:
                currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id = row
                yield seqdata.ReadPair(currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, eval(h1), eval(h2), eval(h3), constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, eval(annotations), fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id, self.analysisfolder)
                #yield seqdata.ReadPair(pairId, header, header, sequence1, sequence2, qual1, qual2,eval(handleCoordinates),clusterId,eval(annotations), fromFastq)

        self.commitAndClose()

    def getReadPairs(self, listOfIds):

        #
        # Imports
        #
        import sys
        from dbs_analysis import seqdata
        import sqlite3
        import time
        import os
          
        inMem = False
        while not inMem:
            try:
                #
                # open connection to database
                #
                self.getConnection()
          
                #
                # faster getting reads from db
                #
                if type(listOfIds[0]) == int: listOfIds = [str(tmp_id) for tmp_id in listOfIds]
                rows = self.c.execute('SELECT id, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id FROM reads WHERE id IN ('+', '.join(listOfIds)+')')#.fetchall()
                for row in rows:
                
                    currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id = row
                    yield seqdata.ReadPair(currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, eval(h1), eval(h2), eval(h3), constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, eval(annotations), fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id, self.analysisfolder)
                
                #
                # alternatively this
                #
                #for readPairId in listOfIds:
                #    #row = self.c.execute('SELECT id,header,sequence1,sequence2,quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq FROM reads WHERE id=?', (int(readPairId), ) ).fetchone()
                #    row = self.c.execute('SELECT id, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos FROM reads WHERE id=?', (int(readPairId), ) ).fetchone()
                #
                #    currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos = row
                #    yield seqdata.ReadPair(currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, eval(h1), eval(h2), eval(h3), constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, eval(annotations), fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos)
                
                self.commitAndClose()
                inMem = True
            except sqlite3.OperationalError as err:
                sys.stderr.write('database pid='+str(os.getpid())+', read loading failed will retry in 1s.\n')
                print err
                time.sleep(1)
 
    def getRuns(self, runTypes):
        
        self.getConnection()
        
        runsInfo = []
        data = self.c.execute('SELECT * FROM runs').fetchall()
        for startTime, command, commandLine, finishedSuccessfully, masterPid in data:
            if command in runTypes: runsInfo.append([startTime, command, commandLine, finishedSuccessfully, masterPid])
        
        self.commitAndClose()
        
        return runsInfo

    @property
    def datadropped(self,):

        #
        # Imports
        #
        import sys
        from dbs_analysis import seqdata
        
        #sys.stderr.write('dropp 1\n')
        self.getConnection()
        #sys.stderr.write('dropp 2\n')
        self.c.execute('SELECT * FROM reads')
        #sys.stderr.write('dropp 3\n')
        self.commitAndClose()
        #sys.stderr.write('dropp 4\n')
        columns = self.c.description
        #sys.stderr.write('dropp 5\n')

        #cursor.execute(query)
        #columns = cursor.description
        #result = []
        #for value in cursor.fetchall():
        #    tmp = {}
        #    for (index,column) in enumerate(value):
        #	tmp[columns[index][0]] = column
        #    result.append(tmp)
        #pprint.pprint(result)

        return bool( len([col[0] for col in columns]) != 32 )

    def dropReadColumns(self,):

        #
        # Imports
        #
        import sys
        from dbs_analysis import seqdata
        
        #
        # open connection to database and drop data
        #
        self.getConnection()
        #self.c.execute("""BEGIN TRANSACTION;")
        self.c.execute("CREATE TEMPORARY TABLE reads_backup(id, header, clusterId, annotations);")
        self.c.execute("INSERT INTO reads_backup SELECT id, header, clusterId, annotations FROM reads;")
        self.c.execute("DROP TABLE reads;")
        self.c.execute("CREATE TABLE reads(id, header, clusterId, annotations);")
        self.c.execute("INSERT INTO reads SELECT id, header, clusterId, annotations FROM reads_backup;")
        self.c.execute("DROP TABLE reads_backup;")
        #self.c.execute("COMMIT;""")
        #self.c.execute("ALTER TABLE reads DROP COLUMN sequenceR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN sequenceR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN qualR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN qualR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN direction")
        #self.c.execute("ALTER TABLE reads DROP COLUMN h1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN h2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN h3")
        #self.c.execute("ALTER TABLE reads DROP COLUMN constructType")
        #self.c.execute("ALTER TABLE reads DROP COLUMN dbsMatch")
        #self.c.execute("ALTER TABLE reads DROP COLUMN dbsSeq")
        #self.c.execute("ALTER TABLE reads DROP COLUMN dbsQual")
        #self.c.execute("ALTER TABLE reads DROP COLUMN mappingFlagR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN refNameR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN refPosR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN mapQR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN cigarR1")
        #self.c.execute("ALTER TABLE reads DROP COLUMN mappingFlagR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN refNameR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN refPosR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN mapQR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN cigarR2")
        #self.c.execute("ALTER TABLE reads DROP COLUMN insertSize")
        #self.c.execute("ALTER TABLE reads DROP COLUMN annotations")
        #self.c.execute("ALTER TABLE reads DROP COLUMN fromFastqId")
        #self.c.execute("ALTER TABLE reads DROP COLUMN r1PositionInFile")
        #self.c.execute("ALTER TABLE reads DROP COLUMN r2PositionInFile")
        #self.c.execute("ALTER TABLE reads DROP COLUMN bamFilePos")
        self.commitAndClose()

    def getAllClustersLoaded(self, analysisfolder, cluster_id_max=False, cluster_id_min=False, skip_singletons=False, skip_htmlTable=False):
        """ function that loads all cluster info available in the database and returns object will all available info, ie there is no need to run cluster.loadInfo() afterwards"""
        
        from dbs_analysis.seqdata import BarcodeCluster
        import sqlite3
        
        success = False
        while not success:
            while self.writeInProgress.value: time.sleep(0.1)
            
            try:
                # get connection to database
                self.getConnection()
                
                # check the table columns
                columnNames = [col[1] for col in self.c.execute('PRAGMA table_info(barcodeClusters)').fetchall()]
                if 'constructTypes' not in columnNames or 'targetInfo' not in columnNames:
                    
                    # get values set clusterinfo and yield cluster
                    if skip_singletons: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations FROM barcodeClusters WHERE clusterTotalReadCount>1')
                    else: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations FROM barcodeClusters')
                    for (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations) in info:
                        constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed = 'None',None,None,None,None,None,'None','None',None,None,False,
                        hetrozygous_positions, high_quality_cluster = None,None
                        if clusterId == None: continue
                        cluster = BarcodeCluster(clusterId,analysisfolder)
                        cluster.setValues(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed, hetrozygous_positions, high_quality_cluster)
                        yield cluster
                
                else: # ie the new columns are present
                    
                    # get values set clusterinfo and yield cluster
                    if cluster_id_min and cluster_id_max:
                        if skip_singletons: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, htmlTable, analyzed, hetrozygous_positions, high_quality_cluster FROM barcodeClusters WHERE (clusterId BETWEEN '+str(cluster_id_min)+' AND '+str(cluster_id_max)+') AND clusterTotalReadCount>1')
                        else: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, htmlTable, analyzed, hetrozygous_positions, high_quality_cluster FROM barcodeClusters WHERE clusterId BETWEEN '+str(cluster_id_min)+' AND '+str(cluster_id_max)+'')
                    else:
                        if skip_singletons:
                            if skip_htmlTable: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, analyzed, hetrozygous_positions, high_quality_cluster FROM barcodeClusters WHERE clusterTotalReadCount>1')
                            else:              info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, htmlTable, analyzed, hetrozygous_positions, high_quality_cluster FROM barcodeClusters WHERE clusterTotalReadCount>1')
                        else: info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, htmlTable, analyzed, hetrozygous_positions, high_quality_cluster FROM barcodeClusters')
                    
                    if not skip_htmlTable:
                        for (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed, hetrozygous_positions, high_quality_cluster) in info:
                            if clusterId == None: continue
                            cluster = BarcodeCluster(clusterId,analysisfolder)
                            cluster.setValues(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed, hetrozygous_positions, high_quality_cluster)
                            yield cluster
                    else:
                        for (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,analyzed, hetrozygous_positions, high_quality_cluster) in info:
                            htmlTable = None
                            if clusterId == None: continue
                            cluster = BarcodeCluster(clusterId,analysisfolder)
                            cluster.setValues(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed, hetrozygous_positions, high_quality_cluster)
                            yield cluster
                self.commitAndClose()
                
                success = True
            except sqlite3.OperationalError: time.sleep(1)

