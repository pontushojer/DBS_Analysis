class Database(object):
    
    def __init__(self, dbPath):
        self.path = dbPath
    
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

	self.c.executemany('INSERT INTO reads VALUES (?,?,?,?,?,?,?,?,?,?)', readsToAdd)
        
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
        readPairs = self.c.execute('SELECT id,header,sequence1,sequence2,quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq FROM reads')
        
	while True:
	    
	    rows = readPairs.fetchmany()#size=readPairs.arraysize)
	    
	    if not rows: break
	    
	    for row in rows:
		pairId,header,sequence1,sequence2,qual1,qual2,handleCoordinates,clusterId,annotations,fromFastq = row
		yield seqdata.ReadPair(pairId, header, header, sequence1, sequence2, qual1, qual2,eval(handleCoordinates),clusterId,eval(annotations), fromFastq)
	
        self.commitAndClose()

    def getReadPairs(self, listOfIds):
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
#        readPairs = self.c.execute('SELECT id,header,sequence1,sequence2,quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq FROM reads WHERE id IN ('+','.join([str(readPairId) for readPairId in listOfIds])+')')
#        
#	while True:
#	    
#	    rows = readPairs.fetchmany()#size=readPairs.arraysize)
#	    
#	    if not rows: break
#	    
#	    for row in rows:
#		pairId,header,sequence1,sequence2,qual1,qual2,handleCoordinates,clusterId,annotations,fromFastq = row
#		yield ReadPair(pairId, header, header, sequence1, sequence2, qual1, qual2,eval(handleCoordinates),clusterId,eval(annotations), fromFastq)
	
	#
	# alternatively this
	#
	for readPairId in listOfIds:
	    row = self.c.execute('SELECT id,header,sequence1,sequence2,quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq FROM reads WHERE id=?', (int(readPairId), ) ).fetchone()
	    pairId,header,sequence1,sequence2,qual1,qual2,handleCoordinates,clusterId,annotations,fromFastq = row
	    yield ReadPair(pairId, header, header, sequence1, sequence2, qual1, qual2,eval(handleCoordinates),clusterId,eval(annotations), fromFastq)
	
        self.commitAndClose()
 
    def getRuns(self, runTypes):
        
        self.getConnection()
        
        runsInfo = []
        data = self.c.execute('SELECT * FROM runs').fetchall()
        for startTime, command, commandLine, finishedSuccessfully, masterPid in data:
            if command in runTypes: runsInfo.append([startTime, command, commandLine, finishedSuccessfully, masterPid])
        
        self.commitAndClose()
        
        return runsInfo

class Results(object,):
    
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
            'constructTypes':None
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
            'constructTypes':'dictionary holding infromation of all the different types of constructs found in the read populations'
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
    
    def __init__(self, analysisfolder):
        """ object holding the settings used for each part of the analysis """
	
        self.analysisfolder = analysisfolder
        import multiprocessing
	import sequences
        
	self.defaultValues = {
	    'debug':False,
	    'uppmaxProject':'b2014005',
	    'parallelProcesses':multiprocessing.cpu_count(),
	    'maxHandleMissMatches':0,
	    'barcodeLength':len(sequences.DBS),
#	    'analysisParts':None,
	    'barcodeMissmatch':0,
#	    'readsPerUmiCutOff':5,
#	    'umiMaxMisMatch':2,
	    'readsPerClusterCutOff':100,
	    'bowtie2Reference':None,
	    'picardPath':None

	}
	self.explenations = {
	    'debug':'Flag for running the scripts in multiprocessing or as single process run [True/False] (default=False)',
	    'uppmaxProject':'Project id used at uppmax for sbatch scripts [bXXXXXXX] (default=b2014005)',
	    'parallelProcesses':'Number of process to run when doing multiprocess parts of analysis (defaul=multiprocessing.cpu_count())',
	    'barcodeLength':'The length of the bead barcode (default='+str(len(sequences.DBS))+')',
#	    'analysisParts':'Parts of the analysis to run specific for each run.',
	    'barcodeMissmatch':'Number of missmatches allowed in the barcode sequence',
	    'maxHandleMissMatches':'Number of missmatches allowed in the handle sequence',
#	    'readsPerUmiCutOff':'Number of reads supporting one UMI for it to passs filters',
#	    'umiMaxMisMatch':'Number of missmatches allowed in the UMI sequence',
	    'readsPerClusterCutOff':'Number of reads supporting a barcode sequence cluster for it to passs filters',
	    'bowtie2Reference':'path to the bowtie 2 reference index',
	    'picardPath':'path to the picard installation to use'
	}
	self.isDefault = {}
	self.setTime = {}

	self.debug = None
	self.uppmaxProject = None
	self.parallelProcesses = None
	self.maxHandleMissMatches = None
	self.barcodeLength = None
#	self.analysisParts = None
	self.barcodeMissmatch = None
#	self.readsPerUmiCutOff = None
#	self.umiMaxMisMatch = None
	self.readsPerClusterCutOff = None
	self.bowtie2Reference = None
	self.picardPath = None
	
	self.setDefaults()

    def setDefaults(self,):
	for variableName, value in self.defaultValues.iteritems():
	    self.__dict__[variableName] = value
	    self.isDefault[variableName] = True
	    self.setTime[variableName] = None
	return 0

    def loadFromDb(self,):
	
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
	import time
	assert variableName in self.explenations,'Error: you are trying to set an undefined variable.\n'
	self.__dict__[variableName]  = value
	self.isDefault[variableName] = False
	self.setTime[variableName]   = time.time()
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
	self.analysisfolder.logfile.write('checking whats in db.\n')
        alreadyInDb = {}
	data = self.analysisfolder.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()
        if data:
            for variableName,default,value,setTime in data:
		self.analysisfolder.logfile.write('processing variable '+variableName+'')
		alreadyInDb[variableName] = True
		
		if variableName in self.__dict__:
		    if default and not self.isDefault[variableName] or setTime < self.setTime[variableName]:
			if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
			self.analysisfolder.logfile.write(', updating from '+str(value)+' to '+str(self.__dict__[variableName])+', old_setTime '+str(setTime)+' new_setTime '+str(self.setTime[variableName])+'.\n')
			self.analysisfolder.database.c.execute('UPDATE settings SET defaultValue=?, value=?, setTime=? WHERE variableName=?', (self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName],variableName))
		    else: self.analysisfolder.logfile.write(' no update needed.\n')
        
        #
        # Add new vars to database
        #
	self.analysisfolder.logfile.write('adding new vars to db:\n')
        for variableName in self.__dict__:
	    if variableName in ['explenations','defaultValues','isDefault','setTime','analysisfolder']:continue
	    if variableName not in alreadyInDb:
		if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
		values = (variableName,self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName])
		self.analysisfolder.database.c.execute('INSERT INTO settings VALUES (?,?,?,?)', values)
		self.analysisfolder.logfile.write('variable '+variableName+' added to db with value '+str(self.__dict__[variableName])+',')
		if self.isDefault[variableName]:self.analysisfolder.logfile.write(' this is the default value.\n')
		else:self.analysisfolder.logfile.write(' non-default value.\n')
	    else: pass#SEAseqPipeLine.logfile.write('variable\t'+variableName+'\talready in db.\n')
        
	self.analysisfolder.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()
        
        return 0

class AnalysisFolder(object):
    """This class represent the analysis outpt folder, hold the structure for it and track the files within"""
    
    def __init__(self, path, logfile=None):

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
        
        import os
	import sqlite3
        if os.path.exists(self.databaseFileName):
            try:
		self.settings.loadFromDb()
		self.results.loadFromDb()
	    except sqlite3.OperationalError:pass

    def create(self, ):

        import os
        for folder in self.folders:
            if not os.path.exists(folder):
                os.mkdir(folder)
        
        self.database.create()
        

    def checkIntegrity(self, ):
        
        import os
        
        if not os.path.isdir(self.path): return 'FAIL: The folder does not exist.'
        
        for folder in self.folders:
            if not os.path.isdir(folder): return 'FAIL: The folder structure is broken.'
        
        return 'PASS'

