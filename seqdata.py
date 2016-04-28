class ReadPair(object):
    
    import misc
    
    def __init__(self, currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile,r2PositionInFile,bamFilePos):

        # original read info
        self.id = currentRead
        self.header   = header
        self.r1Seq      = sequenceR1 	#self.r1Seq      = strip3primN(sequence1)
        self.r2Seq      = sequenceR2 	#self.r2Seq      = strip3primN(sequence2)
        self.r1Qual     = qualR1
        self.r2Qual     = qualR2
        self.fileOrigin = fromFastqId
        self.r1PositionInFile = r1PositionInFile
        self.r2PositionInFile = r2PositionInFile

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
        self.analysisfolder = None
        self.annotations = annotations
        self.overlap = None
        self.brokenSequence = ''
        #self.construct = None
        self.insert= None
        self.clusterId = clusterId

        #graphical representation
        self.outstring = str(self.id)+''.join([' ' for i in range(10-len(str(self.id)))]),
        self.__str__ = self.outstring

        #mapping info
        self.bamFilePos = bamFilePos
        self.mappingFlagR1 = SamFlag(mappingFlagR1)
        self.refNameR1 = refNameR1
        self.refPosR1 = refPosR1
        self.mapQR1 = mapQR1
        self.cigarR1 = cigarR1
        self.mappingFlagR2 = SamFlag(mappingFlagR2)
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
        return (self.id, self.header, None,None,None,None,self.direction,str(self.h1),str(self.h2),str(self.h3),self.construct,self.dbsmatch,self.dbsSeq,self.dbsQual,self.mappingFlagR1.flag,self.refNameR1,self.refPosR1,self.mapQR1,self.cigarR1,self.mappingFlagR2.flag,self.refNameR2,self.refPosR2,self.mapQR2,self.cigarR2,self.insertSize,self.clusterId,str(self.annotations),self.fileOrigin, self.r1PositionInFile,self.r2PositionInFile,self.bamFilePos)

    def fixInsert(self,):
        """ This functions gets the insert sequence depending on what handles where found earlier... I think... """

        # the insert variable has the following layout
                      #r1s  r1q  r2s  r2q
        self.insert = [None,None,None,None]
        # Where
        #r1s is read one sequence
        #r1q is read one quality
        #r2s is read two sequence
        #r2q is read two quality
        
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
        """ this function matches the DBS sequence against the expected pattern """
        
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
        """" function for finding sequenc motifs in the read sequence """
        
        import re
        #matchfunk = hamming_distance

        startPosition = None
        endPosition   = None
        missmatches   = None

        #matchfunk=levenshtein

        #
        # use regular expression to try and find a perfect match
        #
        perfect_match = re.search(matchsequence, readsequence)
        if perfect_match:
            startPosition = perfect_match.start()
            endPosition = perfect_match.end()
            missmatches = 0

        #
        # if perfect match is not found use hamming distance to try and find a "near perfect" match
        #
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
        """ creates the string with colors for each of the reads """
        
        #
        # THIS function needs commenting and revision!!! //EB
        #
        
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
            
            #######################################################################################
            # STRANGE TO RUN THESE FUNCTIONS HERE SHOULD MAYBE BE DONE AT SOME OTHER POINT! //EB  #
            #                                                                                     #
            self.fixInsert()                                                                      #
            #if self.insert: outputSeq+= '--'+self.insert+'--'                                    #
                                                                                                  #
            if self.h3 and self.h1 and self.h2:                                                   #
                outputSeq += ' LOOK HERE!'                                                        #
                self.matchdbs()                                                                   #
                if self.dbsmatch: outputSeq += ' Match!'                                          #
                else: outputSeq += ' No match...'                                                 #
            #outputSeq += '\n'                                                                    #
            #######################################################################################
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
        
        missmatchesAllowed = self.analysisfolder.settings.maxHandleMissMatches

        # look for h1 in read 1
        self.h1 = self.matchSequence(self.r1Seq,sequences.H1,missmatchesAllowed,startOfRead=True)
        startPosition,endPosition,missmatches = self.h1
        if startPosition!=None and startPosition <= 2: self.direction = '1 -> 2'
        if startPosition==None: self.h1 = None
        
        # look for H3 in read one
        if not not self.direction:
            self.h3 = self.matchSequence(self.r1Seq,sequences.H3,missmatchesAllowed,startOfRead=True)
            startPosition,endPosition,missmatches = self.h3
            if startPosition!=None and startPosition <= 2: self.direction = '2 -> 1'
            if startPosition==None: self.h3 = None

        # look for H1 in read 2
        self.h1_in_both_ends = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H1,missmatchesAllowed,startOfRead=True)
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
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H3,missmatchesAllowed,startOfRead=True)
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
            startPosition,endPosition,missmatches = self.matchSequence(checkSeq,revcomp(sequences.H3),missmatchesAllowed)
            if startPosition!=None:
                self.readIntoh3 = True
                self.readIntoh3Coordinates = [startPosition,endPosition,missmatches]
                if not self.h3: self.h3 = True


        if self.direction and not self.h3_in_both_ends:
            
            # find the h2 handle and DBS sequence
            if self.direction == '1 -> 2':
                self.h2 = self.matchSequence(self.r1Seq,revcomp(sequences.H2),missmatchesAllowed)
                if not self.h2[0]: self.h2 = None
                
                if self.h1 and self.h2:
                    self.dbsPrimaryCoordinates = [self.r1Seq,self.h1[1],self.h2[0],self.r1Qual]
                
                if self.h1_in_both_ends: # find secondary h2
                    self.annotations['h2_r2_coordinates'] = self.matchSequence(self.r2Seq,revcomp(H2),missmatchesAllowed)
                    if self.h1in2ndReadCoordinates[0]==0 and self.annotations['h2_r2_coordinates'][0] or (self.h1in2ndReadCoordinates[0] and self.annotations['h2_r2_coordinates'][0]):
                        self.annotations['secondary_dbs_coordinates'] = [self.r2Seq,self.h1in2ndReadCoordinates[1],self.annotations['h2_r2_coordinates'][0]]

            elif self.direction == '2 -> 1':
                self.h2 = self.matchSequence(self.r2Seq,revcomp(sequences.H2),missmatchesAllowed)
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
        """ function that identifies illumina adapter content within the reads and annotates the reads accordingly"""

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

        #
        # setting initial values
        #
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
        import metadata

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

        import random

        self.analysisfolder.database.getConnection()
        clusterIds = list(self.analysisfolder.database.c.execute('SELECT clusterId FROM barcodeClusters').fetchall())
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
        self.annotations = {}

        self.readPairs = []
        self.readPairsById = {}
        self.readPairIdentities = []
        self.readPairsPassFilter = []
        self.readPairsNotPassingFilter = []

        self.contigSequences = []
        self.contigSequencesPassFilter = []
        self.contigSequencesNotPassingFilter = []
        self.nonSingletonContigs = None

        self.analyzed = False

        self.minMAPQ = 0

        self.filesCreated = []

    def loadClusterInfo(self, ):

        import sqlite3, time
        success = False
        while not success:
            while self.analysisfolder.database.writeInProgress.value: time.sleep(0.1)
            try:
                self.analysisfolder.database.getConnection()
                columnNames = [col[1] for col in self.analysisfolder.database.c.execute('PRAGMA table_info(barcodeClusters)').fetchall()]
                if 'constructTypes' not in columnNames:
                    info = self.analysisfolder.database.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations FROM barcodeClusters WHERE clusterId=?',(self.id,)).fetchall()
                    assert len(info) == 1, 'More than one ('+str(len(info))+') clusters found with id '+str(self.id)
                    (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations) = info[0]
                    constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,htmlTable,analyzed = 'None',None,None,None,None,None,'None',None,False
                else:
                    info = self.analysisfolder.database.c.execute('SELECT clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions, htmlTable, analyzed FROM barcodeClusters WHERE clusterId=?',(self.id,)).fetchall()
                    assert len(info) == 1, 'More than one ('+str(len(info))+') clusters found with id '+str(self.id)
                    (clusterId,clusterTotalReadCount,readPairsList,readBarcodeIdentitiesList,clusterBarcodeSequence,clusterBarcodeQuality,contigSequencesList,annotations,constructTypes, readPairsInBamFile, mappedSEReads, SEreadsPassMappingQualityFilter, goodReadPairs, duplicateReadPairs, goodReadPairPositions,htmlTable,analyzed) = info[0]

                self.analysisfolder.database.commitAndClose()

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
                else:
                    self.annotations    = {}

                if constructTypes: self.constructTypes = eval(constructTypes)
                else: self.constructTypes = None
                self.readPairsInBamFile = readPairsInBamFile
                self.mappedSEReads = mappedSEReads
                self.SEreadsPassMappingQualityFilter = SEreadsPassMappingQualityFilter
                self.goodReadPairs = goodReadPairs
                self.duplicateReadPairs = duplicateReadPairs
                if goodReadPairPositions:self.goodReadPairPositions = eval(goodReadPairPositions)
                else:self.goodReadPairPositions=None
                self.tableStr = htmlTable
                self.analyzed = analyzed
                success = True
            except sqlite3.OperationalError: time.sleep(1)

        return 0

    def loadReadPairs(self, ):

        import sqlite3
        import time
        starttime = time.time()
#	  try:
        from misc import Progress
        while self.analysisfolder.database.writeInProgress.value: time.sleep(0.1)
        p = Progress(self.readPairCount, logfile=self.analysisfolder.logfile,unit='cluster_'+str(self.id)+'_reads', mem=True)
        if not self.analysisfolder.database.datadropped:
            for readPair in self.analysisfolder.database.getReadPairs(self.readPairIdsList):
                self.readPairs.append(readPair)
                self.readPairsById[readPair.id] = readPair
                try : p.update()
                except ValueError: pass
#	  except sqlite3.OperationalError: print 'ERROR: BarcodeCluster.loadReadPairs() is giving a sqlite3.OperationalError!!'
        else:
            for readPair in self.analysisfolder.readsdb.getClusterReadPairs(self.id):
                self.readPairs.append(readPair)
                self.readPairsById[readPair.id] = readPair
                try: p.update()
                except ValueError: pass
        print self.id, time.time()-starttime

        return 0

    def generateReadPairsFastq(self, oneFile=False, concatenate=False,r1Part='nobc'):

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
        """ function is not used might be usefull in future though right now only junk """

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

    @property
    def allreadssamechrom(self,):
        tmp = {}
        for readPair in self.readPairs:
            if readPair.mapQR1 >= self.minMAPQ and readPair.mapQR2 >= self.minMAPQ:
                tmp[readPair.refNameR1] = True
                tmp[readPair.refNameR2] = True
        
        count = 0
        for chrom in tmp:
            if chrom != None: count+=1
        
        if count == 1:
            self.chromosome = chrom
            return True
        elif count > 1: return False
        else: return None

    @property
    def allreadssamepos(self,):
        tmp1 = {}
        tmp2 = {}
        for readPair in self.readPairs:
            if readPair.mapQR1 >= self.minMAPQ and readPair.mapQR2 >= self.minMAPQ:
                tmp1[readPair.refPosR1] = True
                tmp2[readPair.refPosR2] = True

        count1 = 0
        for pos in tmp1:
            if pos != None: count1+=1
        count2 = 0
        for pos in tmp2:
            if pos != None: count2+=1

        if count1 == 1 and count2 == 1: return True
        elif count1 > 1 or count2 > 1: return False
        else: return None

    def createBamFile(self,):

        #
        # imports
        #
        from misc import thousandString
        import pysam
        import time
        import subprocess
        import sys
        import os

        #
        # get the reads from the original bam file with all reads and write to new cluster specific bamfile
        #
        bamfile =  pysam.Samfile(self.analysisfolder.dataPath+'/mappedInserts.bam')
        outputBam = pysam.Samfile(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.bam',mode='wb',header=bamfile.header)
        for pair in self.readPairs:
            if pair.bamFilePos:
                bamfile.seek(pair.bamFilePos)
                r1 = bamfile.next()
                r2 = bamfile.next()
                outputBam.write(r1)
                outputBam.write(r2)
        outputBam.close()

        #
        # Sort the bamfile
        #
        picardLogfile = open(self.analysisfolder.temp+'/'+time.strftime("%y%m%d-%H:%M:%S",time.localtime())+'_picardSortSam.cluster_'+str(self.id)+'.log.txt','w')
        self.filesCreated.append(picardLogfile.name)
        command2 = ['java',
                  '-Xmx5g',
                  '-jar',self.analysisfolder.settings.picardPath+'/SortSam.jar',
                  'MAX_RECORDS_IN_RAM=2500000',
                  'INPUT='+self.analysisfolder.temp+'/cluster_'+str(self.id)+'.bam',
                  'OUTPUT='+self.analysisfolder.temp+'/cluster_'+str(self.id)+'.sorted.bam',
                  'SORT_ORDER=coordinate',
                  'CREATE_INDEX=TRUE'
                  ]
        picard = subprocess.Popen(command2,stdout=subprocess.PIPE,stderr=picardLogfile)
        if self.analysisfolder.logfile:
            try: self.analysisfolder.logfile.write('Starting command: '+' '.join(command2)+'\n')
            except ValueError:pass
        errdata = picard.communicate()
        if (picard.returncode != 0 and picard.returncode != None):
            if picard.returncode != 0:
                print '#\n# cmd: '+' '.join( command2 )+'\n#'
                print 'picard view Error code', picard.returncode, errdata, open(picardLogfile.name).read()
                #return 'FAIL'
                sys.exit()
        picardLogfile = open(self.analysisfolder.temp+'/'+time.strftime("%y%m%d-%H:%M:%S",time.localtime())+'_picardMarkDup.cluster_'+str(self.id)+'.log.txt','w')
        self.filesCreated.append(picardLogfile.name)


        #
        # Mark duplicates in the bamfile
        #
        command2 = ['java',
                  '-Xmx5g',
                  '-jar',self.analysisfolder.settings.picardPath+'/MarkDuplicates.jar',
                  'MAX_RECORDS_IN_RAM=2500000',
                  'INPUT='+self.analysisfolder.temp+'/cluster_'+str(self.id)+'.sorted.bam',
                  'OUTPUT='+self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.bam',
                  'METRICS_FILE='+self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.metrics.txt',
                  'CREATE_INDEX=TRUE'
                  ]
        picard = subprocess.Popen(command2,stdout=subprocess.PIPE,stderr=picardLogfile)
        if self.analysisfolder.logfile:
            try:self.analysisfolder.logfile.write('Starting command: '+' '.join(command2)+'\n')
            except ValueError:pass
        errdata = picard.communicate()
        if (picard.returncode != 0 and picard.returncode != None):
            if picard.returncode != 0:
              print '#\n# cmd: '+' '.join( command2 )+'\n#'
              print 'picard view Error code', picard.returncode, errdata, open(picardLogfile.name).read()
              #return 'FAIL'
              sys.exit()
        self.filesCreated.append(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.bam')
        self.filesCreated.append(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.bai')
        self.filesCreated.append(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.metrics.txt')
        
        #
        # remove temporary files
        #
        os.remove(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.bam')
        os.remove(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.sorted.bam')
        os.remove(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.sorted.bai')

        return 0

    def analyze(self):

        #
        # imports
        #
        import time
        from misc import thousandString
        from misc import percentage
        import pysam

        starttime = time.time()

        #
        # Load cluster data and initiate values for counters etc
        #
        print 'Analyzing data in cluster '+str(self.id)
        try: self.analysisfolder.logfile.write('Analyzing data in cluster '+str(self.id)+' ... '+'\n')
        except ValueError:
            #self.analysisfolder.logfile = open(self.analysisfolder.logfile.name,'a')
            #self.analysisfolder.logfile.write('Analyzing data in cluster '+str(self.id)+' ... '+'\n')
            pass

        constructTypes = {}
        readPairsInBamFile = 0
        readPairsInBamFileCheck = 0
        mappedSEReads = 0
        SEreadsPassMappingQualityFilter = 0
        goodReadPairs = 0
        duplicateReadPairs = 0
        goodReadPairsRows = ''
        unmappedReadPairsRows = ''
        duplicateReadPairsRows = ''
        goodReadPairPositions = {}

        #
        # Load the read pairs from database
        #
        try: self.analysisfolder.logfile.write('Loading reads for cluster '+str(self.id)+' ... '+'\n')
        except ValueError: pass
        self.loadReadPairs()

        #
        # build the bamfile with read mappings
        #
        try:self.analysisfolder.logfile.write('Creating bamfiles for cluster_'+str(self.id)+' ... '+'\n')
        except ValueError: pass
        self.createBamFile()
        try:self.analysisfolder.logfile.write('Bamfiles ready for cluster_'+str(self.id)+'.'+'\n')
        except ValueError: pass
        bamfile = pysam.Samfile(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.markedDuplicates.bam')

        #
        # checking constructypes in read population...
        #
        for pair in self.readPairs:
            if pair.bamFilePos: readPairsInBamFile+=1
            try: constructTypes[ pair.construct ] += 1
            except KeyError: constructTypes[ pair.construct ] = 1

        #
        # Parse the bamfile
        #
        try:self.analysisfolder.logfile.write('Making reads table forcluster '+str(self.id)+'.\n')
        except ValueError: pass
        for alignedReadRead in bamfile:

            #
            # get the refernce names for both reads in the read pair
            #
            if not alignedReadRead.is_unmapped: r1ReferenceName = bamfile.getrname(alignedReadRead.tid)
            else: r1ReferenceName = 'NA'
            if not alignedReadRead.mate_is_unmapped: r2ReferenceName = bamfile.getrname(alignedReadRead.rnext)
            else: r2ReferenceName = 'NA'
  
            #
            # Filters
            #
            nicePair = bool( r1ReferenceName == r2ReferenceName and not alignedReadRead.is_unmapped and not alignedReadRead.mate_is_unmapped )
            passMappingQuality = bool(alignedReadRead.mapq >= int(self.analysisfolder.settings.mapqCutOff))
            if not alignedReadRead.is_unmapped: mappedSEReads += 1
            if passMappingQuality: SEreadsPassMappingQualityFilter += 1
  
            #
            # To do paired end checks choose read one
            #
            if alignedReadRead.is_read1: 

                readPairsInBamFileCheck += 1

                #
                # Create a HTM table row for quick vizualisation later - THIS COULD BE CHANGED TO CSV TO ANEABLE SOME NICER D3JS STUFF
                #
                row = '<tr><td>'
                if nicePair and passMappingQuality and not alignedReadRead.is_duplicate:
                    row += '<font color="green">'
                elif nicePair and passMappingQuality:
                    row += '<font color="blue">'
                else:
                    row += '<font color="red">'
                row += alignedReadRead.qname+ '</td><td>'+str(alignedReadRead.flag)+'</td><td>'+str(r1ReferenceName)+'</td><td>'+str(r2ReferenceName)+'</td>'
                if alignedReadRead.pos: row += '<td>'+thousandString(alignedReadRead.pos)+'</td>'
                else:row += '<td>'+str(alignedReadRead.pos)+'</td>'
                if alignedReadRead.pnext: row += '<td>'+thousandString(alignedReadRead.pnext)+'</td>'
                else:row += '<td>'+str(alignedReadRead.pnext)+'</td>'
                row += '<td>'+str(abs(alignedReadRead.isize))+'</td><td>'+str(alignedReadRead.mapq)+'</td><td>'+str(alignedReadRead.cigar)+'</td><td>'+str(alignedReadRead.is_proper_pair)+'</td>'
                row += '</tr>'
                if nicePair and passMappingQuality and not alignedReadRead.is_duplicate:
                    goodReadPairsRows += row
                elif nicePair and passMappingQuality:
                    duplicateReadPairsRows += row
                else:
                    unmappedReadPairsRows += row

                #
                # Save coordinate of reads that mapp in nice pairs
                #
                if nicePair and passMappingQuality:
                    goodReadPairs += 1
                    if alignedReadRead.is_duplicate: duplicateReadPairs += 1
                    else:
                        try:
                            goodReadPairPositions[r1ReferenceName].append(min(alignedReadRead.pos,alignedReadRead.pnext))
                        except KeyError:
                            goodReadPairPositions[r1ReferenceName]   =   [min(alignedReadRead.pos,alignedReadRead.pnext)]

        #
        # assert that the database and bamfile counts match
        #
        assert readPairsInBamFile == readPairsInBamFileCheck, '## ERROR ## : The number of reads in the bamfile is not correct!\n'

        #
        # Make the full html table in one string
        #
        headerRow = '<tr>'+'<th>header</th>'+'<th>flags</th>'+'<th>refchrom R1</th>'+'<th>refchrom R2</th>'+'<th>pos R1</th>'+'<th>pos R2</th>'+'<th>insertsize</th>'+'<th>mapQ</th>'+'<th>CIGAR</th>'+'<th>ProperPair</th>'+'</tr>'
        self.tableStr = '<table>' +headerRow+ goodReadPairsRows + duplicateReadPairsRows + unmappedReadPairsRows + '</table>'

        #
        # save some of the stats for later
        #
        self.constructTypes = constructTypes
        self.readPairsInBamFile = readPairsInBamFile
        self.mappedSEReads = mappedSEReads
        self.SEreadsPassMappingQualityFilter = SEreadsPassMappingQualityFilter
        self.goodReadPairs = goodReadPairs
        self.duplicateReadPairs = duplicateReadPairs
        self.goodReadPairPositions = goodReadPairPositions

        try:self.analysisfolder.logfile.write('Cluster '+str(self.id)+' analyzed in '+str(round(time.time()-starttime,2))+' seconds '+'\n')
        except ValueError: pass
        self.analyzed = True

    def updatedb(self,doUpdate=True,returnTuple=False):
        if doUpdate:
            with self.analysisfolder.database.lock:
                self.analysisfolder.database.writeInProgress.value = True
                import sqlite3
                import time
                updated = False
                while not updated:
                    try: 
                        self.analysisfolder.database.getConnection()
                        self.analysisfolder.database.c.execute('PRAGMA table_info(barcodeClusters)')
                        columnNames = [col[1] for col in self.analysisfolder.database.c.fetchall()]
                        if 'constructTypes' not in columnNames:
                            try: self.analysisfolder.logfile.write('Creating columns in table barcodeClusters in database \n')
                            except ValueError: pass
                            self.analysisfolder.database.getConnection()
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column constructTypes string")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column readPairsInBamFile integer")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column mappedSEReads integer")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column SEreadsPassMappingQualityFilter integer")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column goodReadPairs integer")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column duplicateReadPairs integer")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column goodReadPairPositions string")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column htmlTable string")
                            self.analysisfolder.database.c.execute("alter table barcodeClusters add column analyzed BOOLEAN")
                        self.analysisfolder.database.c.execute('UPDATE barcodeClusters SET annotations=?, constructTypes=?,readPairsInBamFile=?, mappedSEReads=?, SEreadsPassMappingQualityFilter=?, goodReadPairs=?, duplicateReadPairs=?, goodReadPairPositions=?, htmlTable=?, analyzed=? WHERE clusterId=?',(str(self.annotations),str(self.constructTypes),self.readPairsInBamFile,self.mappedSEReads,self.SEreadsPassMappingQualityFilter,self.goodReadPairs,self.duplicateReadPairs,str(self.goodReadPairPositions),self.tableStr,self.analyzed,self.id))
                        self.analysisfolder.database.commitAndClose()
                        self.analysisfolder.database.writeInProgress.value = False
                        updated = True
                        print self.id, updated
                    except sqlite3.OperationalError: time.sleep(1)
        if returnTuple:
            return (str(self.annotations),str(self.constructTypes),self.readPairsInBamFile,self.mappedSEReads,self.SEreadsPassMappingQualityFilter,self.goodReadPairs,self.duplicateReadPairs,str(self.goodReadPairPositions),self.tableStr,self.analyzed,self.id)

    def generateHtmlSummary(self):
        #
        # output for viewer
        #
        from misc import percentage
        from misc import thousandString
        outStr = 'Cluster id='+str(self.id)+', with barcode '+str(self.barcodeSequence)+' has a total of '+str(self.readPairCount)+' read pairs.'+'<br><br>\n'
        outStr += '    ' + str(percentage(self.readPairsInBamFile,self.readPairCount))+'% of the read pairs are in the bamfile ('+str(thousandString(self.readPairsInBamFile))+' of '+str(thousandString(self.readPairCount))+'read pairs ) .<br>\n'
        #outStr += str(self.constructTypes)+'.<br>\n'
        outStr += str(percentage(self.mappedSEReads, self.readPairsInBamFile*2))+'% of the reads (SE) are mapped ('+str(thousandString(self.mappedSEReads))+' of '+str(thousandString(self.readPairsInBamFile*2))+' reads in the bamfile).<br>\n'
        outStr += str(percentage(self.SEreadsPassMappingQualityFilter,self.readPairsInBamFile*2))+'% of the reads (SE) are mapped with mapQ GT q'+   str(self.analysisfolder.settings.mapqCutOff)+'  ('+str(thousandString(self.SEreadsPassMappingQualityFilter))+' of '+str(thousandString(self.readPairsInBamFile*2))+' reads in the bamfile).<br>\n'
        outStr += str(percentage(self.goodReadPairs,self.readPairsInBamFile))+'% of the read Pairs are mapped nice and proper ('+str(thousandString(self.goodReadPairs))+' of '+str(thousandString(self.readPairsInBamFile))+' reads in the bamfile).<br>\n'
        outStr += str(percentage(self.duplicateReadPairs,self.goodReadPairs))+'% of the mapped read pairs are duplicates ('+str(thousandString(self.duplicateReadPairs))+' of '+str(thousandString(self.goodReadPairs))+' mapped reads pairs).<br>\n'
        return outStr

    def makeGoodReadPairPositionsPlots(self,outputfilename = None):
        """ nice and old function currently not used maybe in the future :) """

        import pysam
        chromSizes = {}
        for chrom in pysam.Samfile(self.analysisfolder.dataPath+'/mappedInserts.bam').header['SQ']: chromSizes[chrom['SN']] = chrom['LN']
        
        #
        # Make plos for the chromosomes present
        #
        try: self.analysisfolder.logfile.write('Making plots for cluster '+str(self.id)+'.\n')
        except ValueError: pass
        import matplotlib.pyplot as plt
        import operator
        for chrom, data in self.goodReadPairPositions.iteritems():
  
            binSize = int(int(chromSizes[chrom])/100)
            chrom=str(chrom)
  
            plots = []
            fig, axes = plt.subplots(1, sharex=True)
            plots.append(axes.hist(data,range(0,int(chromSizes[chrom]),binSize),label='chrom '+str(chrom)))#,histtype='step'))
  
            handles, labels = axes.get_legend_handles_labels()
            hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
            handles2, labels2 = zip(*hl)
            axes.legend(handles2, labels2,loc=0,fontsize='small')
  
            axes.set_xlabel('coordinate')
            axes.set_ylabel('# Number of reads')
  
            axes.set_xlim([0,int(chromSizes[chrom])])
  
            #
            # save to file
            #
            if outputfilename:
                plt.savefig(outputfilename,dpi=300,bbox_inches='tight')
                self.filesCreated.append(outputfilename)
            else:
                #plt.savefig(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.histogram.'+chrom+'.pdf',dpi=300,bbox_inches='tight')
                #self.filesCreated.append(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.histogram.'+chrom+'.pdf')
                plt.savefig(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.histogram.'+chrom+'.png',dpi=300,bbox_inches='tight')
                self.filesCreated.append(self.analysisfolder.temp+'/cluster_'+str(self.id)+'.histogram.'+chrom+'.png')

    def removeAllFiles(self):
        import os
        for filename in self.filesCreated:
            if os.path.exists(filename):os.remove(filename)

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

def bits(i,n): 
    return tuple((0,1)[i>>j & 1] for j in xrange(n-1,-1,-1))

class SamFlag():

    def __init__(self, flag):
        """ Function that takes a flag and gets the explanation"""
        ########################-Format Description-#########################
        #	Flag	Chr	Description											#
        #1	0x0001	p	the read is paired in sequencing					#
        #2	0x0002	P	the read is mapped in a proper pair					#
        #3	0x0004	u	the query sequence itself is unmapped				#
        #4	0x0008	U	the mate is unmapped								#
        #5	0x0010	r	strand of the query (1 for reverse)					#
        #6	0x0020	R	strand of the mate									#
        #7	0x0040	1	the read is the first read in a pair				#
        #8	0x0080	2	the read is the second read in a pair				#
        #9	0x0100	s	the alignment is not primary						#
        #10	0x0200	f	the read fails platform/vendor quality checks		#
        #11	0x0400	d	the read is either a PCR or an optical duplicate	#
        #####################################################################
  
        self.readtype = None
        self.properpair = None
        self.mapped = None
        self.matemapped = None
        self.strand = None
        self.mate_strand = None
        self.readnum = None
        self.primaryalignment = None
        self.passfilter = None
        self.pcrduplicate = None

        if flag == None:
            self.flag = flag
            return None
        self.flag = int(flag)
        binary = bits(self.flag,11)
        output = ''
        self.readnum = None
        try:
            if int(binary[-1]):	self.readtype = "PE"
            else:			self.readtype = 'SE'
            if int(binary[-2]):	self.properpair = True
            else:			self.properpair = False
            if int(binary[-3]):	self.mapped = False
            else:			self.mapped = True
            if int(binary[-4]):	self.matemapped = False
            else:			self.matemapped = True
            if int(binary[-5]):	self.strand = 'rev'
            else:			self.strand = 'fwd'
            if int(binary[-6]):	self.mate_strand = "rev"
            else:			self.mate_strand = "fwd"
            if int(binary[-7]):	self.readnum = 1
            else:			pass#self.read1 = False
            if int(binary[-8]):	self.readnum = 2
            else:			pass#self.read2 = False
            if int(binary[-9]):	self.primaryalignment = False;self.output += 'not primary alignment\t'
            else:			self.primaryalignment = True
            if int(binary[-10]):	self.passfilter=False;self.output += 'QC failed\t'
            else:			self.passfilter=True
            if int(binary[-11]):	self.pcrduplicate=True;self.output += 'is PCR duplicate'
            else:			self.pcrduplicate=False;
        except ValueError:
            output += '.'
        return None