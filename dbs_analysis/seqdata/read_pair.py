class ReadPair(object):
    """ object that represent an illumina read pair
    """
    
    import dbs_analysis.misc as misc
    
    def __init__(self, currentRead, header, sequenceR1, sequenceR2, qualR1, qualR2, direction, h1, h2, h3, constructType, dbsMatch, dbsSeq, dbsQual, mappingFlagR1, refNameR1, refPosR1, mapQR1, cigarR1, mappingFlagR2, refNameR2, refPosR2, mapQR2, cigarR2,insertSize, clusterId, annotations, fromFastqId, r1PositionInFile,r2PositionInFile,bamFilePos,individual_id,analysisfolder):

        # original read info
        # comming from the raw fastq files, except id that is a increasing int for each entry in fastq
        # the sequences (r1Seq & r2Seq) and qualities (r1Qual & r2Qual) are never saved to the database (except as None values) and are thus not usable after the handle identification step, they could theoretically be added to database though will take a lot of exttra space
        self.id = currentRead
        self.header   = header
        self.r1Seq      = sequenceR1 	#self.r1Seq      = strip3primN(sequence1)
        self.r2Seq      = sequenceR2 	#self.r2Seq      = strip3primN(sequence2)
        self.r1Qual     = qualR1
        self.r2Qual     = qualR2
        self.fileOrigin = fromFastqId # this is the id in the database of the fastq file that the read pair comes from
        self.r1PositionInFile = r1PositionInFile # this is the positions in the original files where the read is found
        self.r2PositionInFile = r2PositionInFile

        # handle flags and coordinates
        # theese varibles stores handle coordinates for the three main handles
        self.direction = direction
        self.h1 = h1
        self.h2 = h2
        self.h3 = h3
        self.construct = constructType # the contsructtype dict stores information about which handles were found dunring handle identification

        # these variables are used to find the handles connected to the individual id and the id sequence/barcode
        # the handles are only used during the dbs_handle_identifier step and are never saved to the database and are therefore defined as None here
        self.fwd_primer = None
        self.individual_id_primer = None
        self.individual_id = individual_id

        # dbs flags and coordinates
        self.dbs = None
        self.dbsSeq = dbsSeq # this is the actual dbs sequence part of the read pair
        self.dbsQual = dbsQual
        self.dbsmatch = dbsMatch # flag true or false depensing on if the dbs sequence math the expected pattern
        self.dbsPrimaryCoordinates = None

        # other flags and information
        self.analysisfolder = None # this is a holder for defining an analysis folder link
        self.annotations = annotations # a dict with other info that does not really fit anywhere else
        self.overlap = None
        self.brokenSequence = ''
        #self.construct = None
        self.insert= None # this is a list with the part of the sequence that are between the handles only used in dbs_handle_identifier then saved to various fastq files and used for mapping
        self.clusterId = clusterId # the cluster id that the read pair belong to after clustering inititally None

        # graphical representation, not really used anymore, still kept here if needed in the future
        self.outstring = str(self.id)+''.join([' ' for i in range(10-len(str(self.id)))]),
        self.__str__ = self.outstring

        from dbs_analysis.seqdata import SamFlag
        # mapping info from the bamfile look at pysam and sam/bam file defenition to get info about each field/variable
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
        
        # connection to analysisfolder object
        self.analysisfolder = analysisfolder

    @property
    def databaseTuple(self, ):
        """ returning a tuple that is formated to be added to sql database"""
        # data base has following info:
        #       (id,          header,  sequence1, sequence2, quality1,quality2,handleCoordinates,clusterId,annotation,fromFastq)

        # Dumping the quality and sequence values
        #return  (self.id,     self.r1Header,self.r1Seq,self.r2Seq,self.r1Qual,self.r2Qual,   str(self.handleCoordinates),self.dbs,     str(self.annotations),self.fileOrigin)
        #return  (self.id,     self.r1Header,None,None,None,None,   str(self.handleCoordinates),self.dbs,     str(self.annotations),self.fileOrigin)
        return (self.id, self.header, None,None,None,None,self.direction,str(self.h1),str(self.h2),str(self.h3),self.construct,self.dbsmatch,self.dbsSeq,self.dbsQual,self.mappingFlagR1.flag,self.refNameR1,self.refPosR1,self.mapQR1,self.cigarR1,self.mappingFlagR2.flag,self.refNameR2,self.refPosR2,self.mapQR2,self.cigarR2,self.insertSize,self.clusterId,str(self.annotations),self.fileOrigin, self.r1PositionInFile,self.r2PositionInFile,self.bamFilePos, self.individual_id)

    def fixInsert(self,):
        """ This functions gets the insert sequence depending on what handles where found earlier"""

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
            if self.individual_id and self.fwd_primer:
                self.insert[0] = self.r1Seq[self.fwd_primer[0]:]
                self.insert[2] = self.r1Qual[self.fwd_primer[0]:]
        
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
        
        
        from dbs_analysis.seqdata import UIPAC2REGEXP
        
        if self.analysisfolder.settings.type == 'HLA':   from dbs_analysis.sequences import HLA_DBS as DBS
        elif self.analysisfolder.settings.type == 'WFA': from dbs_analysis.sequences import WFA_DBS as DBS

        if self.dbsPrimaryCoordinates:
            self.dbs = self.dbsPrimaryCoordinates[0][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]
            self.dbsQual = self.dbsPrimaryCoordinates[3][self.dbsPrimaryCoordinates[1]:self.dbsPrimaryCoordinates[2]]

        if self.dbs:

            dbsSeq=self.dbs.replace('.','')
            dbsSeq = dbsSeq.split('  ')
            if len(dbsSeq)!=1:dbsSeq = ''
            else: dbsSeq = dbsSeq[0]

            dbsRegex = UIPAC2REGEXP(DBS)
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
            
    def matchSequence(self, readsequence, matchsequence, maxDistance, matchfunk=misc.hamming_distance, startOfRead=False,breakAtFirstMatch=False):
        """" function for finding sequenc motifs in the read sequence returns a list with [startcoordinate, endcoordinate, number of missmatches]"""
        
        #
        # Imports
        #
        import re
        from dbs_analysis import misc
        # matchfunk = hamming_distance        
        try:
            from dbs_analysis.hamming_cython_solution import hamming_loop
            matchfunk = hamming_loop
            from dbs_analysis import hamming_cython_solution
        except ImportError: pass
        
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
            if startOfRead: goto = 10
            for i in range(goto):
                
                if i+len(matchsequence) <= len(readsequence): dist = matchfunk(matchsequence,readsequence[i:i+len(matchsequence)])
                else: dist = 10001
                # if breakAtFirstMatch and (matchfunk == misc.hamming_distance or matchfunk == hamming_cython_solution.hamming_loop):
                #     if dist >= maxDistance:
                #         return [i,i+len(matchsequence),dist]
                
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

    def makeColoredOut(self, terminalcolor=True, htmlcolor=False):
        """ creates the string with colors for each of the reads """
        
        #
        # THIS function needs commenting and revision!!! //EB
        #
        from dbs_analysis import misc
        from dbs_analysis.seqdata import revcomp
        
        if terminalcolor:
            color = misc.TerminalColors()
            newline = '\n'
        elif htmlcolor:
            color = misc.HtmlColors()
            newline = '<br>'
        
        outputSeq = ''
        if True:
            if self.direction == '1 -> 2' or self.direction == '? -> ?' or self.direction == None: quals = [self.r1Qual,self.r2Qual[::-1]]
            if self.direction == '2 -> 1': quals = [self.r2Qual,self.r1Qual[::-1]]
            outputSeq += 'QualBar'
            for i in quals[0]:
                qNum = ord(i)-33
                if qNum >= 28:
                    outputSeq += color.GreenIntense
                elif qNum >= 20:
                    outputSeq += color.YellowIntense
                elif qNum >= 0:
                    outputSeq += color.RedIntense
                outputSeq += '_' + color.Color_Off
            outputSeq += ' '
            for i in quals[1]:
                qNum = ord(i)-33
                if qNum >= 28:
                    outputSeq += color.GreenIntense
                elif qNum >= 20:
                    outputSeq += color.YellowIntense
                elif qNum >= 0:
                    outputSeq += color.RedIntense
                outputSeq += '_' + color.Color_Off
        outputSeq += newline
        
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
                    outputSeq += self.r1Seq[lastWritten:self.h1[0]]+color.BlueIntense+self.r1Seq[self.h1[0]:self.h1[1]]+color.Color_Off
                    lastWritten = self.h1[1]
                
                if self.h2:
                    outputSeq += self.r1Seq[lastWritten:self.h2[0]]+color.PurpleIntense+self.r1Seq[self.h2[0]:self.h2[1]]+color.Color_Off
                    lastWritten = self.h2[1]
                    
                outputSeq += self.r1Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h3:
                    if  self.h3 == True and self.readIntoh3 == True: outputSeq += '###### SOMETHING HERE #####'
                    else:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h3[1]]+color.YellowIntense+revcomp(self.r2Seq)[len(self.r2Seq)-self.h3[1]:len(self.r2Seq)-self.h3[0]]+color.Color_Off
                        lastWritten = len(self.r2Seq)-self.h3[0]
                    if self.h1_in_both_ends and self.readIntoh3 == None:
                        import sys
                        outputSeq += revcomp(self.r2Seq)[lastWritten:]
                        print ' ########  wowowow! funky buissyness!',outputSeq,str(self.h1in2ndReadCoordinates),'  ##################'
                        sys.exit()
                    
                if self.h1_in_both_ends:
                    
                    if self.annotations['h2_r2_coordinates'][0]:
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]]+color.PurpleIntense+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]]+color.Color_Off
                        lastWritten = len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]
                    
                    if self.h1in2ndReadCoordinates[0] == 0 or (self.h1in2ndReadCoordinates[0] != None and self.h1in2ndReadCoordinates[0] != False):
                        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h1in2ndReadCoordinates[1]]+color.BlueIntense+revcomp(self.r2Seq)[len(self.r2Seq)-self.h1in2ndReadCoordinates[1]:len(self.r2Seq)-self.h1in2ndReadCoordinates[0]]+color.Color_Off
                        lastWritten = len(self.r2Seq)-self.h1in2ndReadCoordinates[0]
                outputSeq += revcomp(self.r2Seq)[lastWritten:]
                    
            elif self.direction == '2 -> 1':
                
                if self.h1:
                    #outputSeq += 'HERE --'+self.r2Seq+' --\n'
                    outputSeq += self.r2Seq[lastWritten:self.h1[0]]+color.BlueIntense+self.r2Seq[self.h1[0]:self.h1[1]]+color.Color_Off
                    lastWritten = self.h1[1]
                
                if self.h2:
                    outputSeq += self.r2Seq[lastWritten:self.h2[0]]+color.PurpleIntense+self.r2Seq[self.h2[0]:self.h2[1]]+color.Color_Off
                    lastWritten = self.h2[1]
                    
                if self.annotations and self.readIntoh3Coordinates != None:
                    outputSeq += self.r2Seq[lastWritten:self.readIntoh3Coordinates[0]]+color.YellowIntense+self.r2Seq[self.readIntoh3Coordinates[0]:self.readIntoh3Coordinates[1]]+color.Color_Off
                    lastWritten = self.readIntoh3Coordinates[1]
                    #self.readIntoh3Coordinates = [startPosition,endPosition,missmatches]
                
                if self.readIntoh3Coordinates == None:
                    outputSeq += self.r2Seq[lastWritten:]
                
                outputSeq += ' '
                lastWritten = 0
                if self.h3 and self.h3 != True:
                    outputSeq += revcomp(self.r1Seq)[lastWritten:len(self.r1Seq)-self.h3[1]]+color.YellowIntense+revcomp(self.r1Seq)[len(self.r1Seq)-self.h3[1]:len(self.r1Seq)-self.h3[0]]+color.Color_Off
                    lastWritten = len(self.r1Seq)-self.h3[0]
                    #if self.h1_in_both_ends:
                    #    import sys
                    #    print 'wowowow!'
                    #    sys.exit()
                    
                #if self.h1_in_both_ends:
                #    
                #    if self.annotations['h2_r2_coordinates'][0]:
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]]+color.PurpleIntense+revcomp(self.r2Seq)[len(self.r2Seq)-self.annotations['h2_r2_coordinates'][1]:len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]]+color.Color_Off
                #        lastWritten = len(self.r2Seq)-self.annotations['h2_r2_coordinates'][0]
                #    
                #    if self.h1in2ndReadCoordinates[0] == 0 or (self.h1in2ndReadCoordinates[0] != None and self.h1in2ndReadCoordinates[0] != False):
                #        outputSeq += revcomp(self.r2Seq)[lastWritten:len(self.r2Seq)-self.h1in2ndReadCoordinates[1]]+color.BlueIntense+revcomp(self.r2Seq)[len(self.r2Seq)-self.h1in2ndReadCoordinates[1]:len(self.r2Seq)-self.h1in2ndReadCoordinates[0]]+color.Color_Off
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
        """ function that uses the matchSequence function to find the handles defined in sequences
        """
        from dbs_analysis.seqdata import revcomp
        from dbs_analysis.seqdata import UIPAC2REGEXP
        if self.analysisfolder.settings.type == 'HLA':
            from dbs_analysis.sequences import HLA_H1 as H1
            from dbs_analysis.sequences import HLA_H2 as H2
            from dbs_analysis.sequences import HLA_H3 as H3
            from dbs_analysis.sequences import HLA_DBS as DBS
        elif self.analysisfolder.settings.type == 'WFA':
            from dbs_analysis.sequences import WFA_H1 as H1
            from dbs_analysis.sequences import WFA_H2 as H2
            from dbs_analysis.sequences import WFA_H3 as H3
            from dbs_analysis.sequences import WFA_DBS as DBS
        from dbs_analysis.sequences import IND_HANDLE_1
        from dbs_analysis.sequences import IND_HANDLE_2

        # set direction to None (ie. not identified)
        self.direction = None
        
        missmatchesAllowed = self.analysisfolder.settings.maxHandleMissMatches

        import re
        perfect_read_regex = re.compile('^'+H1+UIPAC2REGEXP(DBS)+revcomp(H2)+'[AGTCN]+'+revcomp(H3)+'$')
        if perfect_read_regex.match(self.r1Seq+revcomp(self.r2Seq)):
            self.direction  = '1 -> 2'
            self.h1 = [0,len(H1),0]
            self.h2 = [len(H1)+len(DBS),len(H1)+len(DBS)+len(H2),0]
            self.h3 = [0,len(H3),0]
            self.dbsPrimaryCoordinates = [self.r1Seq,self.h1[1],self.h2[0],self.r1Qual]
            self.construct = 'constructOK'
            
            self.h1_in_both_ends = None
            self.h3_in_both_ends = None
            self.h1in2ndReadCoordinates = None
            self.h3in2ndReadCoordinates = None
            self.readIntoh3 = None
            self.readIntoh3Coordinates = None
            
            if self.analysisfolder.settings.IndexReferenceTsv:
                self.individual_id_primer = self.matchSequence(self.r1Seq,IND_HANDLE_1,missmatchesAllowed,breakAtFirstMatch=True)
                self.fwd_primer           = self.matchSequence(self.r1Seq,IND_HANDLE_2,missmatchesAllowed,breakAtFirstMatch=True)
                if self.individual_id_primer[0] and self.fwd_primer[0]: self.individual_id = self.r1Seq[self.individual_id_primer[1]:self.fwd_primer[0]]
                else: self.individual_id = None
            else:
                self.individual_id_primer = None
                self.fwd_primer           = None
                self.individual_id = None
            return ''


        # look for h1 in read 1
        self.h1 = self.matchSequence(self.r1Seq,H1,missmatchesAllowed,startOfRead=True,breakAtFirstMatch=True)
        startPosition,endPosition,missmatches = self.h1
        if startPosition!=None and startPosition <= 2: self.direction = '1 -> 2'
        if startPosition==None: self.h1 = None
        
        # look for H3 in read one
        if not self.direction:
            self.h3 = self.matchSequence(self.r1Seq,H3,missmatchesAllowed,startOfRead=True,breakAtFirstMatch=True)
            startPosition,endPosition,missmatches = self.h3
            if startPosition!=None and startPosition <= 2: self.direction = '2 -> 1'
            if startPosition==None: self.h3 = None

        # look for H3 in read two
        self.h3_in_both_ends = None
        startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,H3,missmatchesAllowed,startOfRead=True,breakAtFirstMatch=True)
        if startPosition!=None and startPosition <= 2:
            if self.h3:
                self.h3_in_both_ends = True
                self.h3in2ndReadCoordinates = [startPosition,endPosition,missmatches]
            else:
                self.h3_in_both_ends = False
                self.direction = '1 -> 2'
                self.h3 = [startPosition,endPosition,missmatches]

        # look for H1 in read 2
        self.h1_in_both_ends = None
        if not (self.h3 and self.direction == '1 -> 2'):
            startPosition,endPosition,missmatches = self.matchSequence(self.r2Seq,H1,missmatchesAllowed,startOfRead=True,breakAtFirstMatch=True)
            if startPosition!=None and startPosition <= 2:
                if self.h1:
                    self.h1_in_both_ends = True
                    self.h1in2ndReadCoordinates = [startPosition,endPosition,missmatches]
                else:
                    self.direction = '2 -> 1'
                    self.h1 = [startPosition,endPosition,missmatches]
                    self.h1in2ndReadCoordinates = self.h1
                
        # look for readinto h3
        checkSeq = None
        self.readIntoh3 = None
        self.readIntoh3Coordinates = None
        if self.direction == '1 -> 2': checkSeq = self.r1Seq
        elif self.direction == '2 -> 1': checkSeq = self.r2Seq
        if checkSeq:
            startPosition,endPosition,missmatches = self.matchSequence(checkSeq,revcomp(H3),missmatchesAllowed,breakAtFirstMatch=True)
            if startPosition!=None:
                self.readIntoh3 = True
                self.readIntoh3Coordinates = [startPosition,endPosition,missmatches]
                if not self.h3: self.h3 = True

        # look for individual index in read 1 //HA
        self.individual_id = None
        if self.direction and not self.h3_in_both_ends:
            
            # find the h2 handle and DBS sequence
            if self.direction == '1 -> 2':
                self.h2 = self.matchSequence(self.r1Seq,revcomp(H2),missmatchesAllowed,breakAtFirstMatch=True)
                if not self.h2[0]: self.h2 = None
                
                if self.h1 and self.h2:
                    self.dbsPrimaryCoordinates = [self.r1Seq,self.h1[1],self.h2[0],self.r1Qual]
                    
                    # look for individual id
                    self.individual_id_primer = self.matchSequence(self.r1Seq,IND_HANDLE_1,missmatchesAllowed,breakAtFirstMatch=True)
                    self.fwd_primer           = self.matchSequence(self.r1Seq,IND_HANDLE_2,missmatchesAllowed,breakAtFirstMatch=True)
                    if self.individual_id_primer[0] and self.fwd_primer[0]: self.individual_id = self.r1Seq[self.individual_id_primer[1]:self.fwd_primer[0]]
                    else: self.individual_id = None
                
                if self.h1_in_both_ends: # find secondary h2
                    self.annotations['h2_r2_coordinates'] = self.matchSequence(self.r2Seq,revcomp(H2),missmatchesAllowed,breakAtFirstMatch=True)
                    if self.h1in2ndReadCoordinates[0]==0 and self.annotations['h2_r2_coordinates'][0] or (self.h1in2ndReadCoordinates[0] and self.annotations['h2_r2_coordinates'][0]):
                        self.annotations['secondary_dbs_coordinates'] = [self.r2Seq,self.h1in2ndReadCoordinates[1],self.annotations['h2_r2_coordinates'][0]]

            elif self.direction == '2 -> 1':
                self.h2 = self.matchSequence(self.r2Seq,revcomp(H2),missmatchesAllowed,breakAtFirstMatch=True)
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

    def getFromFiles(self,):
        
        """ function for retreaving data from files on disc not the database
        currently only the fastqs maybe bam aswell later on"""
        
        from dbs_analysis.seqdata import revcomp
        from dbs_analysis import metadata
        
        tmp = self.analysisfolder.database.c.execute('SELECT filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength FROM fastqs WHERE filePairId == '+str(self.fileOrigin)).fetchone()
        filePairId,fastq1,fastq2,readCount,addedToReadsTable,minReadLength = tmp
        
        #
        # get r1 from fastq
        #
        if fastq1.split('.')[-1] in ['gz','gzip']: #
            import gzip                            # open thr fastq file
            fastq1 = gzip.open(fastq1)             # 
        else: fastq1 = open(fastq1)                #
        fastq1.seek(self.r1PositionInFile) # go to the read position
        tmp = fastq1.readline().rstrip().split(' ')[0] # get the header
        assert self.header == tmp, 'ERROR: '+str(self.header)+' != '+str(tmp) # assert that it is the correct read
        self.r1Seq = fastq1.readline().rstrip() # get the sequence
        fastq1.readline().rstrip() # trash the "+"
        self.r1Qual = fastq1.readline().rstrip() # get the quality

        #
        # get r2 from fastq
        #
        if fastq2.split('.')[-1] in ['gz','gzip']: #
            import gzip                            # open thr fastq file
            fastq2 = gzip.open(fastq2)             # 
        else: fastq2 = open(fastq2)                #
        fastq2.seek(self.r2PositionInFile) # go to the read position
        tmp = fastq2.readline().rstrip().split(' ')[0] # get the header
        assert self.header == tmp, 'ERROR: '+str(self.header)+' != '+str(tmp) # assert that it is the correct read
        self.r2Seq = fastq2.readline().rstrip() # get the sequence
        fastq2.readline().rstrip() # trash the "+"
        self.r2Qual = fastq2.readline().rstrip() # get the quality

