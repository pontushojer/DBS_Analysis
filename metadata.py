class Database(object):
    
    def __init__(self, dbPath):
        self.path = dbPath

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
        try: self.conn = sqlite3.connect(self.path)
        except sqlite3.OperationalError:
            print 'ERROR: Trouble with the database, plase check your commandline.'
            sys.exit()
        self.c = self.conn.cursor()
    
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
        import misc
        
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
                    if fastq in filePair:
                        message = 'ERROR: '+fastq+' already in the database.\nExiting after error.'
                        print message
                        if logfile: logfile.write(message+'\n')
                        sys.exit(1)
        #
        # if not in the database add a new row
        #
        if logfile: logfile.write('Getting readcount for file'+fastq1+' ... \n')
        readCount = misc.bufcount(fastq1)/4 #one read is four lines
        if logfile: logfile.write('...done. The file has '+str(readCount)+' reads.\n')
        addedToReadsTable = False#SEAseqPipeLine.startTimeStr
        minReadLength = 'NA'

        if filePairIds: filePairId = max(filePairIds)+1
        else: filePairId = 0
        values = (filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength)
        self.c.execute('INSERT INTO fastqs VALUES (?,?,?,?,?,?)', values)
        
        self.commitAndClose()
        
        return 0
   
    def addReads(self, readsToAdd):

        #
        # Imports
        #
        import sys
        
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
        import seqdata
        
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
                yield seqdata.ReadPair(currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, eval(h1), eval(h2), eval(h3), constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, eval(annotations), fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id)
                #yield seqdata.ReadPair(pairId, header, header, sequence1, sequence2, qual1, qual2,eval(handleCoordinates),clusterId,eval(annotations), fromFastq)

        self.commitAndClose()

    def getReadPairs(self, listOfIds):

        #
        # Imports
        #
        import sys
        import seqdata
        import sqlite3
        import time
          
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
                rows = self.c.execute('SELECT id, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id FROM reads WHERE id IN ('+', '.join(listOfIds)+')').fetchall()
                for row in rows:
                
                    currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id = row
                    yield seqdata.ReadPair(currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, eval(h1), eval(h2), eval(h3), constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, eval(annotations), fromFastqId, r1PositionInFile, r2PositionInFile, bamFilePos, individual_id)
                
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
            except sqlite3.OperationalError: time.sleep(1)
 
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
        import seqdata
        
        self.getConnection()
        self.c.execute('SELECT * FROM reads')
        self.commitAndClose()
        columns = self.c.description

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
        import seqdata
        
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

    def getAllClustersLoaded(self, analysisfolder):
        """ function that loads all cluster info available in the database and returns object will all available info, ie there is no need to run cluster.loadInfo() afterwards"""
        
        from seqdata import BarcodeCluster
        
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
                    info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations FROM barcodeClusters')
                    for (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations) in info:
                        constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed = 'None',None,None,None,None,None,'None','None',None,None,False,
                        if clusterId == None: continue
                        cluster = BarcodeCluster(clusterId,analysisfolder)
                        cluster.setValues(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed)
                        yield cluster
                
                else: # ie the new columns are present
                    
                    # get values set clusterinfo and yield cluster
                    info = self.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, targetInfo,individual_ID_dictionary, htmlTable, analyzed FROM barcodeClusters')
                    for (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed) in info:
                        if clusterId == None: continue
                        cluster = BarcodeCluster(clusterId,analysisfolder)
                        cluster.setValues(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,targetInfo,individual_ID_dictionary,htmlTable,analyzed)
                        yield cluster
                
                self.commitAndClose()
                
                success = True
            except sqlite3.OperationalError: time.sleep(1)

class Results(object,):
    
    #
    # this object has the same basic structure as settings see the settings object for comments
    # in the future this might be better handled by inheretance
    #
    
    def __init__(self, analysisfolder):
        """ object holding the results of the analysis """

        self.analysisfolder = analysisfolder
        
        self.defaultValues = {
            'totalReadCount':None,
            'uniqueBarcodeSequences':None,
            'readPairsHasBarcode':None,
            'readPairsAreIlluminaAdapters':None,
            'barcodeClusterCount':None,
            'singeltonBarcodeClusters':None,
            'minR1readLength':None,
            'minR2readLength':None,
            'readsWithDbsPatternMatch':None,
            'constructTypes':None,
            'bt2AlignmentRate':None,
            'alignmentRateQ20':None
        }
        self.explenations = {
            'totalReadCount':'The number of reads totally included in the analysis.',
            'uniqueBarcodeSequences':'The total number of uniqu barcode sequences identified.',
            'readPairsHasBarcode':'Total number of readpairs where a barcode could be identified.',
            'readPairsAreIlluminaAdapters':'Total number of readpairs where a illumina adapter sequence could be identified.',
            'barcodeClusterCount':'Total number of barcode clusters identified.',
            'singeltonBarcodeClusters':'Total number of barcode clusters of only one read pair identified.',
            'minR1readLength':'Minimum read length found in infiles',
            'minR2readLength':'Minimum read length found in infiles',
            'readsWithDbsPatternMatch':'dictrionary holding the counts for the dbs matching statistics of the reapopulation',
            'constructTypes':'dictionary holding infromation of all the different types of constructs found in the read populations',
            'bt2AlignmentRate':'the observed rate of aligned reads',
            'alignmentRateQ20':'rate of SE reads with mappingQ >= 20'
        }
        self.isDefault = {}
        self.setTime = {}
  
        self.totalReadCount = None
        self.uniqueBarcodeSequences = None
        self.readPairsHasBarcode = None
        self.readPairsAreIlluminaAdapters = None
        self.barcodeClusterCount = None
        self.singeltonBarcodeClusters = None
        self.minR1readLength = None
        self.minR2readLength = None
        self.readsWithDbsPatternMatch = None
        self.constructTypes = None
        self.bt2AlignmentRate = None
        self.alignmentRateQ20 = None
  
        self.setDefaults()

    def setDefaults(self,):
        for resultName, value in self.defaultValues.iteritems():
            self.__dict__[resultName] = value
            self.isDefault[resultName] = True
            self.setTime[resultName] = None
        return 0

    def loadFromDb(self,):
        
        #
        # Get the connection
        #
        self.analysisfolder.database.getConnection()
        
        #
        # Select data
        #
        data = self.analysisfolder.database.c.execute('SELECT resultName,defaultValue,value,setTime FROM results').fetchall()
        
        #
        # Parse data and add to object __dict__
        #
        if data:
            for resultName,default,value,setTime in data:
                self.__dict__[resultName]  = value
                self.isDefault[resultName] = default
                self.setTime[resultName]   = setTime
        
        #
        # close connection
        #
        self.analysisfolder.database.commitAndClose()

    def setResult(self,resultName,value):
        import time
        assert resultName in self.explenations,'Error: you are trying to set a value for an undefined result.\n'
        self.__dict__[resultName]  = value
        self.isDefault[resultName] = False
        self.setTime[resultName]   = time.time()
        return 0

    def saveToDb(self,):
        
        #
        # imports
        #
        import time
        
        #
        # get connection
        #
        self.analysisfolder.database.getConnection()
        
        #
        # Look whats already in database, update it if older or default and set what is not
        #
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('checking results in db.\n')
        alreadyInDb = {}
        data = self.analysisfolder.database.c.execute('SELECT resultName,defaultValue,value,setTime FROM results').fetchall()
        if data:
            for resultName,default,value,setTime in data:
                if self.analysisfolder.logfile: self.analysisfolder.logfile.write('processing result '+resultName+'')
                alreadyInDb[resultName] = True
                
                if resultName in self.__dict__:
                    if default and not self.isDefault[resultName] or setTime < self.setTime[resultName]:
                        if type(self.__dict__[resultName]) in [dict,list]: self.__dict__[resultName] = str(self.__dict__[resultName])
                        if self.analysisfolder.logfile: self.analysisfolder.logfile.write(', updating from '+str(value)+' to '+str(self.__dict__[resultName])+', old_setTime '+str(setTime)+' new_setTime '+str(self.setTime[resultName])+'.\n')
                        self.analysisfolder.database.c.execute('UPDATE results SET defaultValue=?, value=?, setTime=? WHERE resultName=?', (self.isDefault[resultName],self.__dict__[resultName],self.setTime[resultName],resultName))
                    else:
                        if self.analysisfolder.logfile: self.analysisfolder.logfile.write(' no update needed.\n')
        
        #
        # Add new vars to database
        #
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('adding new results to db:\n')
        for resultName in self.__dict__:
            if resultName in ['explenations','defaultValues','isDefault','setTime','analysisfolder']:continue
            if resultName not in alreadyInDb and not self.isDefault[resultName]:
                  if type(self.__dict__[resultName]) in [dict,list]: self.__dict__[resultName] = str(self.__dict__[resultName])
                  values = (resultName,self.isDefault[resultName],self.__dict__[resultName],self.setTime[resultName])
                  self.analysisfolder.database.c.execute('INSERT INTO results VALUES (?,?,?,?)', values)
                  if self.analysisfolder.logfile: self.analysisfolder.logfile.write('result '+resultName+' added to db with value '+str(self.__dict__[resultName])+'\n')
                  #if self.isDefault[resultName]:SEAseqPipeLine.logfile.write(' this is the default value.\n')
                  #else:SEAseqPipeLine.logfile.write(' non-default value.\n')
            else: pass#SEAseqPipeLine.logfile.write('variable\t'+resultName+'\talready in db.\n')
        
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()
        
        return 0

class Settings(object,):
    """ The settings objects stores and saves variables in the database in a nice and controlled manner"""
    
    def __init__(self, analysisfolder):
        """ object holding the settings used for each part of the analysis """
        
        #
        # NOTE: All variables needs to be defined in the Defaults,Explenation and Variable sections below!
        #       adding and removing variables can easily be done by modifying these three sections
        #       without need to modify anything else
        #
        
        self.analysisfolder = analysisfolder
        import multiprocessing
        import sequences
        
        #
        # Default values for all settings
        #
        self.defaultValues = {
            'debug':False,
            'uppmaxProject':'b2014005',
            'parallelProcesses':multiprocessing.cpu_count(),
            'maxHandleMissMatches':0,
            'barcodeLength':len(sequences.DBS),
            #'analysisParts':None,
            'barcodeMissmatch':0,
            #'readsPerUmiCutOff':5,
            #'umiMaxMisMatch':2,
            'readsPerClusterCutOff':0,
            'bowtie2Reference':None,
            'picardPath':None,
            'mapqCutOff':0,
            'minPairsPerCluster':2,
            'targetRegionBed':None

        }
        
        #
        # Explenation strings for all settings
        #
        self.explenations = {
            'debug':'Flag for running the scripts in multiprocessing or as single process run [True/False] (default=False)',
            'uppmaxProject':'Project id used at uppmax for sbatch scripts [bXXXXXXX] (default=b2014005)',
            'parallelProcesses':'Number of process to run when doing multiprocess parts of analysis (defaul=multiprocessing.cpu_count())',
            'barcodeLength':'The length of the bead barcode (default='+str(len(sequences.DBS))+')',
            #'analysisParts':'Parts of the analysis to run specific for each run.',
            'barcodeMissmatch':'Number of missmatches allowed in the barcode sequence',
            'maxHandleMissMatches':'Number of missmatches allowed in the handle sequence',
            #'readsPerUmiCutOff':'Number of reads supporting one UMI for it to passs filters',
            #'umiMaxMisMatch':'Number of missmatches allowed in the UMI sequence',
            'readsPerClusterCutOff':'Number of reads supporting a barcode sequence cluster for it to passs filters',
            'bowtie2Reference':'path to the bowtie 2 reference index',
            'picardPath':'path to the picard installation to use',
            'mapqCutOff':'filter all reads with mapping quality less than this (default='+str(self.defaultValues['mapqCutOff'])+')',
            'minPairsPerCluster':'minimum number of read pairs supporting a cluster for it to be included in analysis (default 2)',
            'targetRegionBed':'A bedfile defining regions on the reference that will be used as a targets during analysis such as coverage stats and variant calling (default = None)'
        }
        
        #
        # containers needed to store information about latest time of seting update to deceide if update is needed
        #
        self.isDefault = {}
        self.setTime = {}
        
        #
        # variables defenition so that they exist in self.__dict__
        #
        self.debug = None
        self.uppmaxProject = None
        self.parallelProcesses = None
        self.maxHandleMissMatches = None
        self.barcodeLength = None
        #self.analysisParts = None
        self.barcodeMissmatch = None
        #self.readsPerUmiCutOff = None
        #self.umiMaxMisMatch = None
        self.readsPerClusterCutOff = None
        self.bowtie2Reference = None
        self.picardPath = None
        self.mapqCutOff = None
        self.minPairsPerCluster=None
        self.targetRegionBed = None
        
        # set the default values
        self.setDefaults()

    def setDefaults(self,):
        """ sets the default valeus to the variables using self.__dict__ and matching self.default  """
        for variableName, value in self.defaultValues.iteritems():
            self.__dict__[variableName] = value
            self.isDefault[variableName] = True
            self.setTime[variableName] = None
        return 0

    def loadFromDb(self,):
        """ Load any info from the database settings table and store it this object
        """
      
        #
        # Get the connection
        #
        self.analysisfolder.database.getConnection()

        #
        # Select data
        #
        data = self.analysisfolder.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()

        #
        # Parse data and add to object __dict__
        #
        if data:
            for variableName,default,value,setTime in data:
                self.__dict__[variableName]  = value
                self.isDefault[variableName] = default
                self.setTime[variableName]   = setTime

        #
        # close connection
        #
        self.analysisfolder.database.commitAndClose()

    def setVariable(self,variableName,value):
        """ update the variable value and the set time etc """
        
        import time
        
        # check that the variable name is valid
        assert variableName in self.explenations,'Error: you are trying to set an undefined variable.\n'
        
        # set the values
        self.__dict__[variableName]  = value
        self.isDefault[variableName] = False
        self.setTime[variableName]   = time.time()
        
        return 0

    def saveToDb(self,):
        """ saves the variables in the object to the database table so that they can be loaded in the next executable or process etc...
        """

        #
        # imports
        #
        import time

        #
        # get connection
        #
        self.analysisfolder.database.getConnection()

        #
        # Look whats already in database, update it if older or default and set what is not
        #
        self.analysisfolder.logfile.write('checking whats in db.\n')
        alreadyInDb = {}
        data = self.analysisfolder.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()
        if data:

            for variableName,default,value,setTime in data:
                self.analysisfolder.logfile.write('processing variable '+variableName+'')
                alreadyInDb[variableName] = True

                if variableName in self.__dict__:
                    # check if the variable is the default value or needs update based on time of the set action
                    if default and not self.isDefault[variableName] or setTime < self.setTime[variableName]:
                        
                        # convert any dicts or lists to strings
                        if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
                        
                        # send a message to the log and update database
                        self.analysisfolder.logfile.write(', updating from '+str(value)+' to '+str(self.__dict__[variableName])+', old_setTime '+str(setTime)+' new_setTime '+str(self.setTime[variableName])+'.\n')
                        self.analysisfolder.database.c.execute('UPDATE settings SET defaultValue=?, value=?, setTime=? WHERE variableName=?', (self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName],variableName))
                    
                    else: self.analysisfolder.logfile.write(' no update needed.\n')

        #
        # Add new vars to database
        #
        self.analysisfolder.logfile.write('adding new vars to db:\n')
        for variableName in self.__dict__:
            
            # if not variable but a known container of other types of info such as explenations or settimes
            if variableName in ['explenations','defaultValues','isDefault','setTime','analysisfolder']:continue
            
            # new sets ie not updates which were already handled above
            if variableName not in alreadyInDb:
                
                # to str conversion like earlier
                if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
                
                # update the database
                values = (variableName,self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName])
                self.analysisfolder.database.c.execute('INSERT INTO settings VALUES (?,?,?,?)', values)
                
                # some info to logfile
                self.analysisfolder.logfile.write('variable '+variableName+' added to db with value '+str(self.__dict__[variableName])+',')
                if self.isDefault[variableName]:self.analysisfolder.logfile.write(' this is the default value.\n')
                else:self.analysisfolder.logfile.write(' non-default value.\n')
            
            else: pass#SEAseqPipeLine.logfile.write('variable\t'+variableName+'\talready in db.\n')

        self.analysisfolder.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()
        
        return 0

class AnalysisFolder(object):
    """This class represent the analysis outpt folder, hold the structure for it and track the files within it is also the container that enables other info such as settings to be sent around among functions in a controlled way"""
    
    def __init__(self, path, logfile=None):
        """ initiates the object and sets the paths to all files relative to commandline input
        also creates the objects database and settings which will be used extensively during the analysis
        """

        #
        # imports
        #
        import os
        import sqlite3

        # Folders
        self.path = path
        self.logpath = self.path+'/logfiles'
        self.rawdataPath = self.path+'/rawData'
        self.dataPath = self.path+'/data'# output/input files
        self.temp = self.path+'/temp'# temp
        
        self.folders = [self.path,self.logpath,self.rawdataPath,self.dataPath,self.temp]

        #Files
        self.databaseFileName  = self.path+'/database.db'
        self.report = self.path+'/report'# report
        self.fastq_outfile1 = self.dataPath+'/inserts.r1.fastq'
        self.fastq_outfile2 = self.dataPath+'/inserts.r2.fastq'
        self.fastq_outfile3 = self.dataPath+'/inserts.singlets.fastq'
        self.coloredReadMasking = self.dataPath+'/coloredReadMasking'
        self.dbsfastq = self.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq'
        
        self.filenames = [self.databaseFileName,self.report]
        
        # objects
        self.database  = Database(self.databaseFileName)
        self.settings = Settings(self)
        self.results = Results(self)
        
        # if database already exists load the info from it
        if os.path.exists(self.databaseFileName):
            try:
                self.settings.loadFromDb()
                self.results.loadFromDb()
            except sqlite3.OperationalError: pass

    def create(self, ):
        """ This functions creates the folder-structure and the database """

        import os
        for folder in self.folders:
            if not os.path.exists(folder):
                os.mkdir(folder)
        
        self.database.create()
        
    def checkIntegrity(self, ):
        """ A function that checks if all folders are in place """
        
        import os
        
        if not os.path.isdir(self.path): return 'FAIL: The folder does not exist.'
        
        for folder in self.folders:
            if not os.path.isdir(folder): return 'FAIL: The folder structure is broken.'
        
        return 'PASS'