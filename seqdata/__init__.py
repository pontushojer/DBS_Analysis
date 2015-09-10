class ReadPair(object):
    
    import misc
    
    def __init__(self, currentRead, header1, header2, sequence1, sequence2, qual1, qual2,handleCoordinates,clusterId,annotations, fastq1):

        # original read info
        self.id = currentRead
	self.r1Header   = header1
	self.r2Header   = header2
	self.r1Seq      = sequence1 	#self.r1Seq      = strip3primN(sequence1)
	self.r2Seq      = sequence2 	#self.r2Seq      = strip3primN(sequence2)
	self.r1Qual     = qual1
	self.r2Qual     = qual2
	self.fileOrigin = fastq1
	
        # handle flags and coordinates
        self.handleCoordinates = handleCoordinates
        self.h770 = None
        self.h1691= None
        self.h4328= None

        # dbs flags and coordinates
        self.dbs = None
        self.dbsmatch = None
        self.dbsPrimaryCoordinates = None

        # other flags and information
	self.annotations= annotations
        self.overlap = None
        self.brokenSequence = ''
        self.construct = None
        self.insert= None
        
        #graphical representation
        self.outstring = str(self.id)+''.join([' ' for i in range(10-len(str(self.id)))]),
        self.__str__ = self.outstring

    def fixInsert(self,):
        
        if self.direction == '1 -> 2':
            if self.h1691:
                self.insert = self.r1Seq[self.h1691[1]:]
                if 'readinto_h4328_coordinates' in self.annotations:
                    self.insert = self.r1Seq[self.h1691[1]:self.annotations['readinto_h4328_coordinates'][0]]
        elif self.direction == '2 -> 1':
            if self.h1691:
                self.insert = self.r2Seq[self.h1691[1]:]
                if 'readinto_h4328_coordinates' in self.annotations:
                    self.insert = self.r2Seq[self.h1691[1]:self.annotations['readinto_h4328_coordinates'][0]]

        #if self.insert:
        #    self.insert=self.insert.replace('.','')
        #    self.insert = self.insert.split('  ')
        #    if len(self.insert)>1:
        #        if len(self.insert[0])>len(self.insert[1]):
        #            longInsertPart = self.insert[0]
        #            shortInsertPart = self.insert[1]
        #        else:
        #            longInsertPart = self.insert[1]
        #            shortInsertPart = self.insert[0]
        #            
        #        self.insert = longInsertPart
        #    else: self.insert = self.insert[0]

                #start,end,missmatches = self.matchSequence(longInsertPart,shortInsertPart,int(round(len(shortInsertPart)*0.1))+2)
                
                #print longInsertPart[:start] + '-'+longInsertPart[start:end] + '-' + longInsertPart[end:]
                
            #print self.insert
        return 0

    def matchdbs(self,):

        if self.dbsPrimaryCoordinates:
            self.dbs = self.dbsPrimaryCoordinates[0][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]

        if self.dbs:

            dbsSeq=self.dbs.replace('.','')
            dbsSeq = dbsSeq.split('  ')
            if len(dbsSeq)!=1:dbsSeq = ''
            else: dbsSeq = dbsSeq[0]

            dbs = 'BDHVBDHVBDHVBDHVBDHV'
            dbsRegex = UIPAC2REGEXP(dbs)
            import re
            
            if dbsSeq:
                match = re.match(dbsRegex,self.dbs)
                if match:
                    self.dbsmatch = True
                    #print 'woohooo!',dbsSeq
                else:
                    self.dbsmatch = False
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

    def matchHandles(self,):
        
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
                
                if self.h770:
                    outputSeq += self.r1Seq[lastWritten:self.h770[0]]+'\033[34m'+self.r1Seq[self.h770[0]:self.h770[1]]+'\033[0m'
                    lastWritten = self.h770[1]
                
                if self.h1691:
                    outputSeq += self.r1Seq[lastWritten:self.h1691[0]]+'\033[95m'+self.r1Seq[self.h1691[0]:self.h1691[1]]+'\033[0m'
                    lastWritten = self.h1691[1]
                    
                outputSeq += self.r1Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h4328:
                    if  self.h4328 == True and self.annotations['readinto_h4328'] == True: outputSeq += '###### SOMETHING HERE #####'
                    else:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h4328[1]]+'\033[93m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.h4328[1]:len(self.r2Seq)-self.h4328[0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.h4328[0]
                    if self.annotations['h770_in_both_ends'] and not'readinto_h4328' in self.annotations:
                        import sys
                        outputSeq += revcomp(self.r2Seq)[lastWritten:]
                        print ' ########  wowowow! funky buissyness!',outputSeq,str(self.annotations['h770_r2_coordinates']),'  ##################'
                        sys.exit()
                    
                if self.annotations['h770_in_both_ends']:
                    
                    if self.annotations['h1691_r2_coordinates'][0]:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][1]]+'\033[95m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][0]
                    
                    if self.annotations['h770_r2_coordinates'][0] == 0 or (self.annotations['h770_r2_coordinates'][0] != None and self.annotations['h770_r2_coordinates'][0] != False):
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h770_r2_coordinates'][1]]+'\033[34m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h770_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h770_r2_coordinates'][0]]+'\033[0m'
                        lastWritten = len(self.r2Seq)-self.annotations['h770_r2_coordinates'][0]
                outputSeq += revcomp(self.r2Seq)[lastWritten:]
                    
            elif self.direction == '2 -> 1':
                
                if self.h770:
                    #outputSeq += 'HERE --'+self.r2Seq+' --\n'
                    outputSeq += self.r2Seq[lastWritten:self.h770[0]]+'\033[34m'+self.r2Seq[self.h770[0]:self.h770[1]]+'\033[0m'
                    lastWritten = self.h770[1]
                
                if self.h1691:
                    outputSeq += self.r2Seq[lastWritten:self.h1691[0]]+'\033[95m'+self.r2Seq[self.h1691[0]:self.h1691[1]]+'\033[0m'
                    lastWritten = self.h1691[1]
                    
                if self.annotations and 'readinto_h4328_coordinates' in self.annotations:
                    outputSeq += self.r2Seq[lastWritten:self.annotations['readinto_h4328_coordinates'][0]]+'\033[93m'+self.r2Seq[self.annotations['readinto_h4328_coordinates'][0]:self.annotations['readinto_h4328_coordinates'][1]]+'\033[0m'
                    lastWritten = self.annotations['readinto_h4328_coordinates'][1]
                    #self.annotations['readinto_h4328_coordinates'] = [startPosition,endPosition,missmatches]
                
                if 'readinto_h4328_coordinates' not in self.annotations:
                    outputSeq += self.r2Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h4328 and self.h4328 != True:
                    outputSeq += revcomp(self.r1Seq)[lastWritten:len(self.r1Seq)-self.h4328[1]]+'\033[93m'+revcomp(self.r1Seq)[len(self.r1Seq)-self.h4328[1]:len(self.r1Seq)-self.h4328[0]]+'\033[0m'
                    lastWritten = len(self.r1Seq)-self.h4328[0]
                    #if self.annotations['h770_in_both_ends']:
                    #    import sys
                    #    print 'wowowow!'
                    #    sys.exit()
                    
                #if self.annotations['h770_in_both_ends']:
                #    
                #    if self.annotations['h1691_r2_coordinates'][0]:
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][1]]+'\033[95m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][0]]+'\033[0m'
                #        lastWritten = len(self.r2Seq)-self.annotations['h1691_r2_coordinates'][0]
                #    
                #    if self.annotations['h770_r2_coordinates'][0] == 0 or (self.annotations['h770_r2_coordinates'][0] != None and self.annotations['h770_r2_coordinates'][0] != False):
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h770_r2_coordinates'][1]]+'\033[34m'+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h770_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h770_r2_coordinates'][0]]+'\033[0m'
                #        lastWritten = len(self.r2Seq)-self.annotations['h770_r2_coordinates'][0]
                outputSeq += revcomp(self.r1Seq)[lastWritten:]
                
                # 2->1 END
            self.fixInsert()
            #if self.insert: outputSeq+= '--'+self.insert+'--'

            if self.h4328 and self.h770 and self.h1691:
                outputSeq += ' LOOK HERE!'
                self.matchdbs()
                if self.dbsmatch: outputSeq += ' Match!'
                else: outputSeq += ' No match...'
            #outputSeq += '\n'
            
            if self.h770 and self.h1691 and self.h4328:
                self.construct = 'constructOK'
            else: 
                self.construct =''
                if not self.h770: self.construct += ' h770 '
                if not self.h4328: self.construct += ' h4328 '
                if not self.h1691: self.construct += ' h1691'
                if self.h770 and self.h4328 and not self.h1691 and self.direction == '1 -> 2' and len(self.r1Seq)<60: self.construct = ' h1691-SemiOK'
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
        
        # look for h770 in read 1
        self.h770 = self.matchSequence(self.r1Seq,sequences.H770,4)
        startPosition,endPosition,missmatches = self.h770
        if startPosition!=None and startPosition <= 2:
            self.direction = '1 -> 2'
        if startPosition==None: self.h770 = None
        
        # look for H4328 in read one
        self.h4328 = self.matchSequence(self.r1Seq,sequences.H4328,4)
        startPosition,endPosition,missmatches = self.h4328
        if startPosition!=None and not self.direction and startPosition <= 2:
            self.direction = '2 -> 1'
        if startPosition==None: self.h4328 = None

        # look for H770 in read 2
        self.annotations['h770_in_both_ends'] = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H770,4)
        if startPosition!=None and startPosition <= 2:
            if self.h770:
                self.annotations['h770_in_both_ends'] = True
                self.annotations['h770_r2_coordinates'] = [startPosition,endPosition,missmatches]
            else:
                self.direction = '2 -> 1'
                self.h770 = [startPosition,endPosition,missmatches]
                self.annotations['h770_r2_coordinates'] = self.h770
        
        # look for H4328 in read two
        self.annotations['h4328_in_both_ends'] = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,sequences.H4328,4)
        if startPosition!=None and startPosition <= 2:
            if self.h4328:
                self.annotations['h4328_in_both_ends'] = True
                self.annotations['h4328_r2_coordinates'] = [startPosition,endPosition,missmatches]
            else:
                self.annotations['h4328_in_both_ends'] = False
                self.direction = '1 -> 2'
                self.h4328 = [startPosition,endPosition,missmatches]
        
        # look for readinto h4328
        checkSeq = None
        if self.direction == '1 -> 2': checkSeq = self.r1Seq
        elif self.direction == '2 -> 1': checkSeq = self.r2Seq
        if checkSeq:
            startPosition,endPosition,missmatches = self.matchSequence(checkSeq,revcomp(sequences.H4328),4)
            if startPosition!=None:
                self.annotations['readinto_h4328'] = True
                self.annotations['readinto_h4328_coordinates'] = [startPosition,endPosition,missmatches]
                if not self.h4328: self.h4328 = True

        if self.direction and not self.annotations['h4328_in_both_ends']:
            
            # find the h1691 handle and DBS sequence
            if self.direction == '1 -> 2':
                self.h1691 = self.matchSequence(self.r1Seq,revcomp(sequences.H1691),4)
                if not self.h1691[0]: self.h1691 = None
                
                if self.h770 and self.h1691:
                    self.dbsPrimaryCoordinates = [self.r1Seq,self.h770[1],self.h1691[0]]
                
                if self.annotations['h770_in_both_ends']: # find secondary h1691
                    self.annotations['h1691_r2_coordinates'] = self.matchSequence(self.r2Seq,revcomp(H1691),4)
                    if self.annotations['h770_r2_coordinates'][0]==0 and self.annotations['h1691_r2_coordinates'][0] or (self.annotations['h770_r2_coordinates'][0] and self.annotations['h1691_r2_coordinates'][0]):
                        self.annotations['secondary_dbs_coordinates'] = [self.r2Seq,self.annotations['h770_r2_coordinates'][1],self.annotations['h1691_r2_coordinates'][0]]

            elif self.direction == '2 -> 1':
                self.h1691 = self.matchSequence(self.r2Seq,revcomp(sequences.H1691),4)
                if not self.h1691[0]: self.h1691 = None

                if self.h770 and self.h1691:
                    self.dbsPrimaryCoordinates = [self.r2Seq,self.h770[1],self.h1691[0]]

            else:pass
    
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

