class ReadPair(object):
    
    import misc
    
    def __init__(self, currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId):
	
	# original read info
        self.id = currentRead
	self.header   = header
	self.r1Seq      = sequenceR1 	#self.r1Seq      = strip3primN(sequence1)
	self.r2Seq      = sequenceR2 	#self.r2Seq      = strip3primN(sequence2)
	self.r1Qual     = qualR1
	self.r2Qual     = qualR2
	self.fileOrigin = fromFastqId
	
        # handle flags and coordinates
        self.direction = direction
	self.h1 = h1
        self.h2 = h2
        self.h3 = h3
	self.construct = constructType

        # dbs flags and coordinates
        self.dbs = None
	self.dbsSeq = dbsSeq
	self.dbsQual = dbsQual
        self.dbsmatch = dbsMatch
        self.dbsPrimaryCoordinates = None

        # other flags and information
	self.annotations = annotations
        self.overlap = None
        self.brokenSequence = ''
        self.construct = None
        self.insert= None
	self.clusterId = clusterId
        
        #graphical representation
        self.outstring = str(self.id)+''.join([' ' for i in range(10-len(str(self.id)))]),
        self.__str__ = self.outstring
	
	#mapping info
	self.mappingFlagR1 = mappingFlagR1
	self.refNameR1 = refNameR1
	self.refPosR1 = refPosR1
	self.mapQR1 = mapQR1
	self.cigarR1 = cigarR1
	self.mappingFlagR2 = mappingFlagR2
	self.refNameR2 = refNameR2
	self.refPosR2 = refPosR2
	self.mapQR2 = mapQR2
	self.cigarR2 = cigarR2
	self.insertSize = insertSize

    @property
    def databaseTuple(self, ):
	""" returning a tuple that is formated to be added to sql database"""
	# data base has following info:
	#       (id,          header,  sequence1, sequence2, quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq)
	
	# Dumping the quality and sequence values
	#return  (self.id,     self.r1Header,self.r1Seq,self.r2Seq,self.r1Qual,self.r2Qual,   str(self.handleCoordinates),self.dbs,     str(self.annotations),self.fileOrigin)
	#return  (self.id,     self.r1Header,None,None,None,None,   str(self.handleCoordinates),self.dbs,     str(self.annotations),self.fileOrigin)
	return (self.id, self.header, None,None,None,None,self.direction,str(self.h1),str(self.h2),str(self.h3),self.construct,self.dbsmatch,self.dbsSeq,self.dbsQual,self.mappingFlagR1,self.refNameR1,self.refPosR1,self.mapQR1,self.cigarR1,self.mappingFlagR2,self.refNameR2,self.refPosR2,self.mapQR2,self.cigarR2,self.insertSize,self.clusterId,str(self.annotations),self.fileOrigin)

    def fixInsert(self,):
        
        self.insert = [None,None,None,None]
        
        if self.direction == '1 -> 2':
            if self.h2:
                self.insert[0] = self.r1Seq[self.h2[1]:]
                self.insert[2] = self.r1Qual[self.h2[1]:]
                if self.readIntoh3Coordinates != None:
                    self.insert[0] = self.r1Seq[self.h2[1]:self.readIntoh3Coordinates[0]]
                    self.insert[2] = self.r1Qual[self.h2[1]:self.readIntoh3Coordinates[0]]
            if self.h3 and self.h3 != True:
                self.insert[1] = self.r2Seq[self.h3[1]:]
                self.insert[3] = self.r2Qual[self.h3[1]:]
        
        elif self.direction == '2 -> 1':
            if self.h2:
                self.insert[0] = self.r2Seq[self.h2[1]:]
                self.insert[2] = self.r2Qual[self.h2[1]:]
                if self.readIntoh3Coordinates != None:
                    self.insert[0] = self.r2Seq[self.h2[1]:self.readIntoh3Coordinates[0]]
                    self.insert[2] = self.r2Qual[self.h2[1]:self.readIntoh3Coordinates[0]]
            if self.h3 and self.h3 != True:
                self.insert[1] = self.r1Seq[self.h3[1]:]
                self.insert[3] = self.r1Qual[self.h3[1]:]
        
        if self.insert == [None,None,None,None]: self.insert = None
        return 0

    def matchdbs(self,):
        
        import sequences

        if self.dbsPrimaryCoordinates:
            self.dbs = self.dbsPrimaryCoordinates[0][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]
	    self.dbsQual = self.dbsPrimaryCoordinates[3][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]

        if self.dbs:

            dbsSeq=self.dbs.replace('.','')
            dbsSeq = dbsSeq.split('  ')
            if len(dbsSeq)!=1:dbsSeq = ''
            else: dbsSeq = dbsSeq[0]

            dbsRegex = UIPAC2REGEXP(sequences.DBS)
            import re
            
            if dbsSeq:
                match = re.match('^'+dbsRegex+'$',self.dbs)
                if match:
                    self.dbsmatch = True
		    self.dbsSeq = self.dbs
                    #print 'woohooo!',dbsSeq
                else:
                    self.dbsmatch = False
		    self.dbsSeq = False
                    #print 'ohhhnooooo!',dbsSeq
            else:pass#print 'BADSEQUENCE!',dbsSeq
            
    def matchSequence(self, readsequence, matchsequence, maxDistance, matchfunk=misc.hamming_distance, startOfRead=False):
	
	import re
	#matchfunk = hamming_distance

	startPosition = None
	endPosition   = None
	missmatches   = None

        #matchfunk=levenshtein

	perfect_match = re.search(matchsequence, readsequence)
	if perfect_match:
	    startPosition = perfect_match.start()
	    endPosition = perfect_match.end()
	    missmatches = 0
	
	elif int(maxDistance):
	    mindist = [10000,-1]
            goto = len(readsequence)
            if startOfRead: goto = 30
	    for i in range(goto):
		
		if i+len(matchsequence) <= len(readsequence): dist = matchfunk(matchsequence,readsequence[i:i+len(matchsequence)])
		else: dist = 10001
		
		if dist < mindist[0]: mindist =[dist,i]

	    if mindist[0] < int(maxDistance)+1:
		startPosition = mindist[1]
		endPosition = mindist[1]+len(matchsequence)
		missmatches = mindist[0]
	    else:
		startPosition = None
		endPosition = None
		missmatches = None

	return [startPosition,endPosition,missmatches]

    def makeColoredOut(self,):
        
        outputSeq = ''
        if True:
            if self.direction == '1 -> 2' or self.direction == '? -> ?' or self.direction == None: quals = [self.r1Qual,self.r2Qual[::-1]]
            if self.direction == '2 -> 1': quals = [self.r2Qual,self.r1Qual[::-1]]
            outputSeq += 'QualBar'
            for i in quals[0]:
                qNum = ord(i)-33
                if qNum >= 28:
                    outputSeq += '\033[92m'
                elif qNum >= 20:
                    outputSeq += '\033[93m'
                elif qNum >= 0:
                    outputSeq += '\033[91m'
                outputSeq += '_' + '\033[0m'
            outputSeq += ' '
            for i in quals[1]:
                qNum = ord(i)-33
                if qNum >= 28:
                    outputSeq += '\033[92m'
                elif qNum >= 20:
                    outputSeq += '\033[93m'
                elif qNum >= 0:
                    outputSeq += '\033[91m'
                outputSeq += '_' + '\033[0m'
        outputSeq += '\n'
        
        if   self.direction == '1 -> 2':outputSeq += self.direction+' '
        elif self.direction == '2 -> 1':outputSeq += self.direction+' '
        else:
            #outputSeq = '     ? -> ?'+''.join([' ' for i in range(10-len('? -> ?'))]) +' '+ self.r1Seq+' '+'.'.join(['' for i in range(160-len(self.r2Seq)-len(self.r1Seq))])+' '+revcomp(self.r2Seq)
            outputSeq += '? -> ? '+ self.r1Seq+' '+revcomp(self.r2Seq)
            #print outputSeq
        
        if self.direction:
            
            lastWritten = 0

            if self.direction == '1 -> 2':
                
                if self.h1:
                    outputSeq += self.r1Seq[lastWritten:self.h1[0]]+'\033[34m'+self.r1Seq[self.h1[0]:self.h1[1]]+'\033[0m'
                    lastWritten = self.h1[1]
                
                if self.h2:
                    outputSeq += self.r1Seq[lastWritten:self.h2[0]]+'\033[95m'+self.r1Seq[self.h2[0]:self.h2[1]]+'\033[0m'
                    lastWritten = self.h2[1]
                    
                outputSeq += self.r1Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h3:
                    if  self.h3 == True and self.readIntoh3 == True: outputSeq += '###### SOMETHING HERE #####'
                    else:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h3[1]]+'\033[93m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.h3[1]:len(self.r2Seq)-self.h3[0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.h3[0]
                    if self.h1_in_both_ends and self.readIntoh3 == None:
                        import sys
                        outputSeq += revcomp(self.r2Seq)[lastWritten:]
                        print ' ########  wowowow! funky buissyness!',outputSeq,str(self.h1in2ndReadCoordinates),'  ##################'
                        sys.exit()
                    
                if self.h1_in_both_ends:
                    
                    if self.annotations['h2_r2_coordinates'][0]:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]]+'\033[95m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]
                    
                    if self.h1in2ndReadCoordinates[0] == 0 or (self.h1in2ndReadCoordinates[0] != None and self.h1in2ndReadCoordinates[0] != False):
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h1in2ndReadCoordinates[1]]+'\033[34m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.h1in2ndReadCoordinates[1]:len(self.r2Seq)-self.h1in2ndReadCoordinates[0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.h1in2ndReadCoordinates[0]
                outputSeq += revcomp(self.r2Seq)[lastWritten:]
                    
            elif self.direction == '2 -> 1':
                
                if self.h1:
                    #outputSeq += 'HERE --'+self.r2Seq+' --\n'
                    outputSeq += self.r2Seq[lastWritten:self.h1[0]]+'\033[34m'+self.r2Seq[self.h1[0]:self.h1[1]]+'\033[0m'
                    lastWritten = self.h1[1]
                
                if self.h2:
                    outputSeq += self.r2Seq[lastWritten:self.h2[0]]+'\033[95m'+self.r2Seq[self.h2[0]:self.h2[1]]+'\033[0m'
                    lastWritten = self.h2[1]
                    
                if self.annotations and self.readIntoh3Coordinates != None:
                    outputSeq += self.r2Seq[lastWritten:self.readIntoh3Coordinates[0]]+'\033[93m'+self.r2Seq[self.readIntoh3Coordinates[0]:self.readIntoh3Coordinates[1]]+'\033[0m'
                    lastWritten = self.readIntoh3Coordinates[1]
                    #self.readIntoh3Coordinates = [startPosition,endPosition,missmatches]
                
                if self.readIntoh3Coordinates == None:
                    outputSeq += self.r2Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h3 and self.h3 != True:
                    outputSeq += revcomp(self.r1Seq)[lastWritten:len(self.r1Seq)-self.h3[1]]+'\033[93m'+revcomp(self.r1Seq)[len(self.r1Seq)-self.h3[1]:len(self.r1Seq)-self.h3[0]]+'\033[0m'
                    lastWritten = len(self.r1Seq)-self.h3[0]
                    #if self.h1_in_both_ends:
                    #    import sys
                    #    print 'wowowow!'
                    #    sys.exit()
                    
                #if self.h1_in_both_ends:
                #    
                #    if self.annotations['h2_r2_coordinates'][0]:
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]]+'\033[95m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]]+'\033[0m'
                #        lastWritten = len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]
                #    
                #    if self.h1in2ndReadCoordinates[0] == 0 or (self.h1in2ndReadCoordinates[0] != None and self.h1in2ndReadCoordinates[0] != False):
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h1in2ndReadCoordinates[1]]+'\033[34m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.h1in2ndReadCoordinates[1]:len(self.r2Seq)-self.h1in2ndReadCoordinates[0]]+'\033[0m'
                #        lastWritten = len(self.r2Seq)-self.h1in2ndReadCoordinates[0]
                outputSeq += revcomp(self.r1Seq)[lastWritten:]
                
                # 2->1 END
            self.fixInsert()
            #if self.insert: outputSeq+= '--'+self.insert+'--'

            if self.h3 and self.h1 and self.h2:
                outputSeq += ' LOOK HERE!'
                self.matchdbs()
                if self.dbsmatch: outputSeq += ' Match!'
                else: outputSeq += ' No match...'
            #outputSeq += '\n'
            
            #print self.construct
            
            #if self.dbsPrimaryCoordinates:
            #    print 'PRIMARYDBS',
            #    print self.dbsPrimaryCoordinates[0][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]
            #if 'secondary_dbs_coordinates' in self.annotations:
            #    print 'SCNDARYDBS',
            #    print self.annotations['secondary_dbs_coordinates'][0][self.annotations['secondary_dbs_coordinates'][1]:self.annotations['secondary_dbs_coordinates'][2]]

        #if self.overlap: print 'OLP'+outputSeq
        #else:print ''+outputSeq
        self.outputSeq = outputSeq

    def identifyDirection(self,):
        
        import sequences
        
        # set direction to None (ie. not identified)
        self.direction = None
        
        # look for h1 in read 1
        self.h1 = self.matchSequence(self.r1Seq,sequences.H1,4,startOfRead=True)
        startPosition,endPosition,missmatches = self.h1
        if startPosition!=None and startPosition <= 2: self.direction = '1 -> 2'
        if startPosition==None: self.h1 = None
        
        # look for H3 in read one
        if not not self.direction:
            self.h3 = self.matchSequence(self.r1Seq,sequences.H3,4,startOfRead=True)
            startPosition,endPosition,missmatches = self.h3
            if startPosition!=None and startPosition <= 2: self.direction = '2 -> 1'
            if startPosition==None: self.h3 = None

        # look for H1 in read 2
        self.h1_in_both_ends = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H1,4,startOfRead=True)
        if startPosition!=None and startPosition <= 2:
            if self.h1:
                self.h1_in_both_ends = True
                self.h1in2ndReadCoordinates = [startPosition,endPosition,missmatches]
            else:
                self.direction = '2 -> 1'
                self.h1 = [startPosition,endPosition,missmatches]
                self.h1in2ndReadCoordinates = self.h1
        
        # look for H3 in read two
        self.h3_in_both_ends = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H3,4,startOfRead=True)
        if startPosition!=None and startPosition <= 2:
            if self.h3:
                self.h3_in_both_ends = True
                self.h3in2ndReadCoordinates = [startPosition,endPosition,missmatches]
            else:
                self.h3_in_both_ends = False
                self.direction = '1 -> 2'
                self.h3 = [startPosition,endPosition,missmatches]
        
        # look for readinto h3
        checkSeq = None
        self.readIntoh3 = None
        self.readIntoh3Coordinates = None
        if self.direction == '1 -> 2': checkSeq = self.r1Seq
        elif self.direction == '2 -> 1': checkSeq = self.r2Seq
        if checkSeq:
            startPosition,endPosition,missmatches = self.matchSequence(checkSeq,revcomp(sequences.H3),4)
            if startPosition!=None:
                self.readIntoh3 = True
                self.readIntoh3Coordinates = [startPosition,endPosition,missmatches]
                if not self.h3: self.h3 = True


        if self.direction and not self.h3_in_both_ends:
            
            # find the h2 handle and DBS sequence
            if self.direction == '1 -> 2':
                self.h2 = self.matchSequence(self.r1Seq,revcomp(sequences.H2),4)
                if not self.h2[0]: self.h2 = None
                
                if self.h1 and self.h2:
                    self.dbsPrimaryCoordinates = [self.r1Seq,self.h1[1],self.h2[0],self.r1Qual]
                
                if self.h1_in_both_ends: # find secondary h2
                    self.annotations['h2_r2_coordinates'] = self.matchSequence(self.r2Seq,revcomp(H2),4)
                    if self.h1in2ndReadCoordinates[0]==0 and self.annotations['h2_r2_coordinates'][0] or (self.h1in2ndReadCoordinates[0] and self.annotations['h2_r2_coordinates'][0]):
                        self.annotations['secondary_dbs_coordinates'] = [self.r2Seq,self.h1in2ndReadCoordinates[1],self.annotations['h2_r2_coordinates'][0]]

            elif self.direction == '2 -> 1':
                self.h2 = self.matchSequence(self.r2Seq,revcomp(sequences.H2),4)
                if not self.h2[0]: self.h2 = None

                if self.h1 and self.h2:
                    self.dbsPrimaryCoordinates = [self.r2Seq,self.h1[1],self.h2[0],self.r2Qual]

            else:pass
    
        #classify construct type
        if self.h1 and self.h2 and self.h3:
            self.construct = 'constructOK'
        else: 
            self.construct ='missing:'
            if not self.h1: self.construct += ' h1 '
            if not self.h3: self.construct += ' h3 '
            if not self.h2: self.construct += ' h2'
            if self.h1 and self.h3 and not self.h2 and self.direction == '1 -> 2' and len(self.r1Seq)<60: self.construct = ' h2-SemiOK'
 
    def isIlluminaAdapter(self, ):
	
	import math
	
	for i in [15]:#range(15):
	    handleStartPosition,handleEndPosition,handleMissMatches = self.matchSequence(self.r1Seq, 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'[:35-i], int(math.ceil(float(35-i)*0.1)))
	    if handleStartPosition:
		self.annotations['Read1IsIlluminaAdapter'] = str(handleMissMatches)+':'+ str(handleStartPosition)
		break
	
	for i in [15]:#range(15):
	    handleStartPosition,handleEndPosition,handleMissMatches = self.matchSequence(self.r2Seq, 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'[:34-i], int(math.ceil(float(34-i)*0.1)))
	    if handleStartPosition:
		self.annotations['Read2IsIlluminaAdapter'] = str(handleMissMatches)+':'+ str(handleStartPosition)
		break

class BarcodeClusterer(object):
    
    def __init__(self, analysisfolder):
        
        self.analysisfolder = analysisfolder
        self.logfile = self.analysisfolder.logfile

    def generateBarcodeFastq(self, ):

	#
	# imports
	#
	import operator
        import misc
	
	if self.logfile: self.logfile.write('Generating barcode fastq ...\n')
	
	uniqBarcodeSequences = {}
	temporaryDict = {}
	barcodeCounter = 0
	totalReadPairCounter = 0
	qualities = {}
	readPairHasBarcodeCounter = 0
	
	if self.logfile: self.logfile.write('Loading read pairs ...\n')
	progress = misc.Progress(self.analysisfolder.results.totalReadCount, logfile=self.logfile, unit='reads-loaded-from-db', mem=True)
	with progress:
	    for pair in self.analysisfolder.database.getAllReadPairs():
		if pair.dbsSeq:
		    readPairHasBarcodeCounter += 1
		    barcodeSequence = pair.dbsSeq
		    qualities[pair.id] = pair.dbsQual
		    try:            uniqBarcodeSequences[barcodeSequence].append(pair.id)
		    except KeyError:uniqBarcodeSequences[barcodeSequence] = [pair.id]
		progress.update()
	if self.logfile: self.logfile.write('Done.\n')
	self.analysisfolder.results.setResult('uniqueBarcodeSequences',len(uniqBarcodeSequences))
	if self.logfile: self.logfile.write(str(self.analysisfolder.results.uniqueBarcodeSequences)+' uniq barcode sequences found within the read pair population.\n')

	if self.logfile: self.logfile.write('Sorting the barcodes by number of reads/sequence.\n')
	if self.logfile: self.logfile.write('Building the sorting dictionary ...\n')
	for barcode, idList in uniqBarcodeSequences.iteritems():
		try:		temporaryDict[len(idList)].append(barcode)
		except KeyError:temporaryDict[len(idList)] = [barcode]
	
	
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

	#
	# imports
	#
	import subprocess
        import time

	clusteringProgramsLogfile = open(self.analysisfolder.logpath+'/'+time.strftime("%y%m%d-%H:%M:%S",time.localtime())+'_cdHitBarcodeClustering.log.txt','w')
	
	# cluster raw barcodes
	command = ['cd-hit-454',
		    '-i',self.analysisfolder.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq',
		    '-o',self.analysisfolder.dataPath+'/clusteredBarcodeSequences',
		    '-g','1',      # mode
		    '-n','3',      # wrodsize
		    '-M','0',      # memory limit
		    '-t',str(self.analysisfolder.settings.parallelProcesses),      # threads
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
	command = ['cdhit-cluster-consensus',
		self.analysisfolder.dataPath+'/clusteredBarcodeSequences.clstr',
		self.analysisfolder.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq',
		self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus',
		self.analysisfolder.dataPath+'/clusteredBarcodeSequences.aligned'
		]
	if self.logfile: self.logfile.write('Starting command: '+' '.join(command)+'\n')
	ccc = subprocess.Popen(command,stdout=clusteringProgramsLogfile,stderr=subprocess.PIPE )
	errdata = ccc.communicate()
	if ccc.returncode != 0:
		print 'cmd: '+' '.join( command )
		print 'cdhit-cluster-consensus view Error code', ccc.returncode, errdata
		sys.exit()
	clusteringProgramsLogfile.close()
	
	return 0

    def parseBarcodeClusteringOutput(self, readPairsHasBarcode):
	
	import misc
        
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
	consensusFile = open(self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus.fastq')
	clstrFile = open(self.analysisfolder.dataPath+'/clusteredBarcodeSequences.clstr')
	
	#
	# load cluster ids and consensus sequences
	#
	if self.logfile: self.logfile.write('\nLoading barcode clusters and a barcode consesnsus sequences for each cluster ...\n')
	while True:
	    header = consensusFile.readline().rstrip()
	    barcodeSequence = consensusFile.readline().rstrip()
	    junk = consensusFile.readline().rstrip()
	    barcodeQuality = consensusFile.readline().rstrip()
	    if header == '': break
	    totalClusterCount += 1
	    header = header.split('_cluster_')
	    clusterId = int(header[1].split(' ')[0])
	    if header[0][:2] == '@s':
		singletonClusters[clusterId] = {'clusterReadCount':1,'readPairs':[],'identities':[],'clusterBarcodeSequence':barcodeSequence,'clusterBarcodeQuality':barcodeQuality}
		barcodeClusters[clusterId] = singletonClusters[clusterId]
	    elif header[0][:2] == '@c':
		nonSingletonClusters[clusterId] = {'clusterReadCount':int(header[1].split(' ')[2]),'readPairs':[],'identities':[],'clusterBarcodeSequence':barcodeSequence,'clusterBarcodeQuality':barcodeQuality}
		barcodeClusters[clusterId] = nonSingletonClusters[clusterId]
	    else: raise ValueError
	self.analysisfolder.results.setResult('barcodeClusterCount',totalClusterCount)
	self.analysisfolder.results.setResult('singeltonBarcodeClusters',len(singletonClusters))
	if self.logfile: self.logfile.write('A total of '+str(totalClusterCount)+' clusters of barcode sequences were loaded into memory.\n')

	#
	# Load what readpairs are in each cluster
	#
	if self.logfile: self.logfile.write('\nLoading read pair to barcode cluster connections ...\n')
	progress = misc.Progress(readPairsHasBarcode, logfile=self.logfile, unit='reads-loaded', mem=True)
	with progress:
	    for line in clstrFile:
		line = line.rstrip()
		if line[0] == '>':
		    clusterId = int(line.split(' ')[1])
		    continue
		elif line[0] == '0':
		    readId = line.split('>')[1].split('.')[0]
		    identity = 'seed'
		    assert line.split(' ')[-1] == '*', 'Error in file format of clstr file'
		else:
		    readId = line.split('>')[1].split('.')[0]
		    identity = float(line.split(' ')[-1].split('/')[-1].split('%')[0])
		barcodeClusters[clusterId]['readPairs'].append(readId)
		barcodeClusters[clusterId]['identities'].append(identity)
		progress.update()
	if self.logfile: self.logfile.write('All read pair to barcode cluster connections loaded.\n')
	
	return barcodeClusters

    def addBarcodeClusterInfoToDatabase(self, barcodeClusters):

	import misc
        
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
	if self.logfile: self.logfile.write('Converting the data ... \n')
	for clusterId,data in barcodeClusters.iteritems():
	    for readPairId in data['readPairs']: tmpUpdateValues[readPairId]=clusterId
	    addValues.append( (clusterId,data['clusterReadCount'],str(data['readPairs']),str(data['identities']),data['clusterBarcodeSequence'],data['clusterBarcodeQuality'],None,None) )
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
		self.analysisfolder.dataPath+'/clusteredBarcodeSequences.consensus.fastq',
		self.analysisfolder.dataPath+'/clusteredBarcodeSequences.aligned',
		]
	processes = [] 
	if self.logfile: temp = self.logfile
        else: temp = sys.stdout
        for fileName in a:
	    if self.logfile: self.logfile.write('Starting process:'+' gzip -v9 '+fileName+'\n')
	    p=subprocess.Popen(['gzip','-v9',fileName],stdout=temp)
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
	
	import random
	
	self.analysisfolder.database.getConnection()
	clusterIds = self.analysisfolder.database.c.execute('SELECT clusterId FROM barcodeClusters').fetchall()
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
		if clusterId[0] == None: continue
		clusterId = int(clusterId[0])
		cluster = BarcodeCluster(clusterId,self.analysisfolder)
		cluster.loadClusterInfo()
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

class BarcodeCluster(object,):
    
    def __init__(self, clusterId,analysisfolder):
	
	self.id = clusterId
	self.analysisfolder = analysisfolder
	self.readPairCount = None
	self.contigCount = None
	self.barcodeSequence = None
	self.barcodeQuality = None
	self.readPairIdsList = []
	self.contigIdsList = []
	self.annotations = None

	self.readPairs = []
	self.readPairsById = {}
	self.readPairIdentities = []
	self.readPairsPassFilter = []
	self.readPairsNotPassingFilter = []

	self.contigSequences = []
	self.contigSequencesPassFilter = []
	self.contigSequencesNotPassingFilter = []
	self.nonSingletonContigs = None
	
	self.filesCreated = []

    def loadClusterInfo(self, ):
	
	self.analysisfolder.database.getConnection()
	info = self.analysisfolder.database.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations FROM barcodeClusters WHERE clusterId=?',(self.id,)).fetchall()
	self.analysisfolder.database.commitAndClose()
	assert len(info) == 1, 'More than one ('+str(len(info))+') clusters found with id '+str(self.id)
	(clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations) = info[0]
	assert clusterId == self.id
	self.readPairCount      = int(clusterTotalReadCount)
	self.readPairIdsList    = eval(readPairsList)
	self.readPairIdentities = eval(readBarcodeIdentitiesList)
	self.barcodeSequence    = clusterBarcodeSequence
	self.barcodeQuality     = clusterBarcodeQuality
	if contigSequencesList:
	    self.contigIdsList  = eval(contigSequencesList)
	if annotations:
	    self.annotations    = eval(annotations)
	return 0

    def loadReadPairs(self, ):
	
	import sqlite3
#	try:
	for readPair in self.analysisfolder.database.getReadPairs(self.readPairIdsList):
	    self.readPairs.append(readPair)
	    self.readPairsById[readPair.id] = readPair
#	except sqlite3.OperationalError: print 'ERROR: BarcodeCluster.loadReadPairs() is giving a sqlite3.OperationalError!!'
	
	return 0

    def generateReadPairsFastq(self, oneFile=True, concatenate=True,r1Part='nobc'):
	
	#
	# Create files
	#
	if oneFile:
	    fastqFile = open(SEAseqPipeLine.tempFileFolder+'/cluster_'+str(self.id)+'.fastq','w')
	    self.filesCreated.append(fastqFile.name)
	else:
	    fastqFileR1 = open(SEAseqPipeLine.tempFileFolder+'/cluster_'+str(self.id)+'.R1.fastq','w')
	    fastqFileR2 = open(SEAseqPipeLine.tempFileFolder+'/cluster_'+str(self.id)+'.R2.fastq','w')
	    self.filesCreated.append(fastqFileR1.name)
	    self.filesCreated.append(fastqFileR2.name)

	if not self.readPairs: self.loadReadPairs()
	
	#
	# check readcount
	#
	assert len(self.readPairs) == self.readPairCount, 'BarcodeCluster.readcount does not match number of reads for cluster='+self.id
	
	#
	# create output and write to file(s)
	#
	for readPair in self.readPairs:

	    # remove barcode sequence and handle if requested
	    if r1Part == 'full':
		r1Seq  = readPair.r1Seq
		r1Qual = readPair.r1Qual
	    elif r1Part == 'nobc':
		r1Seq  = readPair.r1Seq[ readPair.handleCoordinates[1]:SEAseqPipeLine.results.minR1readLength]
		r1Qual = readPair.r1Qual[readPair.handleCoordinates[1]:SEAseqPipeLine.results.minR1readLength]
	    
	    # concatenated output or not
	    if concatenate:
		r1Header = '@'+str(readPair.id)#readPair.r1Header+' '+str(readPair.annotations)
		for annotation in readPair.annotations: r1Header += ' '+annotation+'='+str(readPair.annotations[annotation])
		r1 = r1Header+'\n'+r1Seq+'CCCCCCCCCCGGGGGGGGGG'+readPair.r2Seq[:SEAseqPipeLine.results.minR2readLength]+'\n'+'+'+'\n'+r1Qual+'CCCCCCCCCCGGGGGGGGGG'+readPair.r2Qual[:SEAseqPipeLine.results.minR2readLength]+'\n'
		r2 = ''
	    else:
		r1 = readPair.r1Header.split(' ')[0]+'/1 '+' '.join(readPair.r1Header.split(' ')[1:])+'\n'+r1Seq+'\n'+'+'+'\n'+r1Qual+'\n'
		r2 = readPair.r2Header.split(' ')[0]+'/2 '+' '.join(readPair.r1Header.split(' ')[1:])+'\n'+readPair.r2Seq[:SEAseqPipeLine.results.minR2readLength]+'\n'+'+'+'\n'+readPair.r2Qual[:SEAseqPipeLine.results.minR2readLength]+'\n'
	    
	    # write to file(s)
	    if oneFile:
		    fastqFile.write(r1+r2)
	    else:
		fastqFileR1.write(r1)
		fastqFileR2.write(r2)
	    
	#
	# close file connection(s)
	#
	if oneFile:
	    fastqFile.close()
	else:
	    fastqFileR1.close()
	    fastqFileR2.close()
	
	return 0

    def saveContigInfoToDatabase(self, ):

	import sqlite3
	try:
	    addValues = []
	    SEAseqPipeLine.database.getConnection()
	    if self.nonSingletonContigs:
		for contigId, data in self.nonSingletonContigs.iteritems():
		    for readPairId in data['readPairs']:
			for annotation, value in self.readPairsById[int(readPairId)].annotations.iteritems(): data['annotations'][annotation]=value
		    #                  contigId, contigTotalReadCount,      readPairsList,     readPairIdentitiesList,    consensusSequence,       consensusQuality,    clusterId,annotations
		    addValues.append( (contigId,data['contigReadCount'],str(data['readPairs']),str(data['identities']),data['consensusSequence'],data['consensusQuality'],self.id,str(data['annotations'])) )
		    self.contigIdsList.append(contigId)
		SEAseqPipeLine.database.c.executemany('INSERT INTO contigs VALUES (?,?,?,?,?,?,?,?)', addValues)
		SEAseqPipeLine.database.c.execute('UPDATE barcodeClusters SET contigSequencesList=? WHERE clusterId=?', (str(self.contigIdsList), self.id))
		#print self.id, self.contigIdsList
		SEAseqPipeLine.database.commitAndClose()
		#print 'Write OK!!'
	except sqlite3.OperationalError: print 'ERROR: BarcodeCluster.saveContigInfoToDatabase() is giving a sqlite3.OperationalError!!'

	return 0

def readPairGenerator(fastq1,fastq2):

    #
    # imports
    #
    import gzip
    import sys
    import misc

    #
    # Loop through infiles
    #
    currentRead = 0
    totalReadcount = 0
    readcount = misc.bufcount(fastq1)/4
    readfromdiskProgress = misc.Progress(readcount, unit='reads-read-from-disk', mem=True)
    with readfromdiskProgress:
        #for filePairId, readcount, fastq1, fastq2 in self.infiles:
        sys.stderr.write(str(currentRead)+' read pairs read from infiles, now starting to read from '+fastq1+'.\n')
        totalReadcount += readcount
        
        #
        # Open the files
        #
        if fastq1.split('.')[-1] in ['gz','gzip']: file1 = gzip.open(fastq1)
        else: file1 = open(fastq1,'r')
        if fastq2.split('.')[-1] in ['gz','gzip']: file2 = gzip.open(fastq2)
        else: file2 = open(fastq2,'r')
        
        while 'NOT EOFError':
            try:
                header1 = file1.readline().rstrip()
                header2 = file2.readline().rstrip()
                sequence1 = file1.readline().rstrip()
                sequence2 = file2.readline().rstrip()
                trash = file1.readline().rstrip()
                trash = file2.readline().rstrip()
                qual1 = file1.readline().rstrip()
                qual2 = file2.readline().rstrip()
                if not header1: break
                currentRead += 1
                # data base has following info:
                #    (id,header,sequence1,sequence2,quality1,quality2,barcodeSequence,clusterId,annotation,fromFastq)
                barcodeSequence = None
                clusterId = None
                annotations = {}
                if len(sequence1) == 0:
                    sequence1 = 'NoSequence'
                    qual1 =  'NoSequence'
                if len(sequence2) == 0:
                    sequence2 = 'NoSequence'
                    qual2 =  'NoSequence'
                pair = ReadPair(currentRead, header1, header2, sequence1, sequence2, qual1, qual2,barcodeSequence,clusterId,annotations, 0)#fastq1)
                readfromdiskProgress.update()
                yield pair#self.currentRead, header1, header2, sequence1, sequence2, qual1, qual2, fastq1
            except EOFError: break
            #if currentRead >= 10000:break # to testrun on subset
        #assert totalReadcount == currentRead, 'Error while reading infiles: Read count after file '+fastq1+' is '+str(currentRead)+' should theoretically be '+str(totalReadcount)+'.\n'
        sys.stderr.write('Reached the end of '+fastq1+'.\n')
    grandTotal = totalReadcount
    sys.stderr.write(str(grandTotal)+' read pairs read from infiles.\n')
    #SEAseqPipeLine.results.setResult('totalReadCount',grandTotal)

def revcomp(string):
	''' Takes a sequence and reversecomplements it'''
	complementary = comp(string)
	return complementary[::-1]

def comp(string):
	''' Takes a sequence and complements it'''
	complement = {'A':'T','T':'A',
				  'C':'G','G':'C',
				  'N':'N',
				  'R':'Y','Y':'R',
				  'K':'M','M':'K',
				  'B':'V','V':'B',
				  'D':'H','H':'D',
				  }
	compseq = "".join([complement.get(nt.upper(), '') for nt in string])
	return compseq

def strip3primN(string):
    
    while string and string[-1] == 'N': string = string[:-1]
    
    return string

def uipac(bases, back='uipac'): #U	Uracil NOT SUPPORTED!!!
	if back == 'uipac':
		if 'N' in bases: return 'N'
		uniqbases={}
		for i in bases:
			uniqbases[i]=True
		bases = uniqbases.keys()
		if 'U' in bases: sys.stderr.write('WARNING in function "uipac(bases)": Uracil NOT SUPPORTED!')
		if len(bases)==1:
			if 'A' in bases: return 'A' #A	Adenine
			if 'C' in bases: return 'C' #C	Cytosine
			if 'G' in bases: return 'G' #G	Guanine
			if 'T' in bases: return 'T' #T	Thymine
			#U	Uracil NOT SUPPORTED!!!
		if len(bases)==2:
			if 'A' in bases and 'G' in bases: return 'R' #R	Purine (A or G)
			if 'C' in bases and 'T' in bases: return 'Y' #Y	Pyrimidine (C, T, or U)
			if 'A' in bases and 'C' in bases: return 'M' #M	C or A
			if 'T' in bases and 'G' in bases: return 'K' #K	T, U, or G
			if 'A' in bases and 'T' in bases: return 'W' #W	T, U, or A
			if 'C' in bases and 'G' in bases: return 'S' #S	C or G
		if len(bases)==3:
			if 'C' in bases and 'T' in bases and 'G' in bases: return 'B' #B	C, T, U, or G (not A)
			if 'A' in bases and 'T' in bases and 'G' in bases: return 'D' #D	A, T, U, or G (not C)
			if 'A' in bases and 'T' in bases and 'C' in bases: return 'H' #H	A, T, U, or C (not G)
			if 'A' in bases and 'C' in bases and 'G' in bases: return 'V' #V	A, C, or G (not T, not U)
		if len(bases)==4:
			return 'N' #N	Any base (A, C, G, T, or U)
	elif back == 'bases':
		if   bases == 'R': return ['A','G'] 
		elif bases == 'Y': return ['C','T']
		elif bases == 'M': return ['A','C']
		elif bases == 'K': return ['G','T']
		elif bases == 'W': return ['A','T']
		elif bases == 'S': return ['C','G']
		elif bases == 'B': return ['C','T','G']
		elif bases == 'D': return ['A','T','G']
		elif bases == 'V': return ['A','C','G']
		elif bases == 'H': return ['A','C','T']
		elif bases == 'N': return ['A','G','T','C']

def UIPAC2REGEXP(string):
    return string.replace('R','[AG]').replace('Y','[CT]').replace('S','[GC]').replace('W','[AT]').replace('K','[GT]').replace('M','[AC]').replace('B','[CGT]').replace('D','[AGT]').replace('H','[ACT]').replace('V','[ACG]').replace('N','[AGTC]')

