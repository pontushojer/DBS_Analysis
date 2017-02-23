
import itertools

#
# CREATE ALL BASECOMBOS FOR THE REPEAT CHECK
#
bases = [i for i in "AGTC"]
twoBaseCombos   = [ ''.join(combo) for combo in itertools.product('AGTC',repeat=2)]
threeBaseCombos = [ ''.join(combo) for combo in itertools.product('AGTC',repeat=3)]

class Handle():
    """ Handle object
    """

    def __init__(self,
            handle_length = 30,
            handle_upper_gc = 0.70,
            handle_lower_gc = 0.40,
            handle_min_tm = 68,
            handle_max_tm = 100,
            kmerLength = 5,
            minHD = 3
        ):
      
        import random
        self.sequence = ''.join([random.choice("AGCT") for i in range(handle_length)])
        self.logstring = ''
        self.dead = False
        self.resonFordeath = None
        self.upper_gc = handle_upper_gc
        self.lower_gc = handle_lower_gc
        self.min_tm = handle_min_tm
        self.max_tm = handle_max_tm
        self.kmerLength = kmerLength
        self.minHD = minHD

    @property
    def GCcontent(self):
       return round( (self.sequence.count('G')+self.sequence.count('C')) / float(len(self.sequence)),2)

    @property
    def basicTm(self):
       return self.sequence.count('A')*2+self.sequence.count('T')*2+self.sequence.count('G')*4+self.sequence.count('C')*4

    @property
    def betterTm(self):
      return 64.9 + 41*(self.sequence.count('G')+self.sequence.count('C')-16.4)/len(self.sequence)
   
    @property
    def isRepeat(self):
        """ Function identifying repeated patterns in handle sequence
        """
        import re
         
        for repeatUnit in bases:
            # look for repeated bases in stretches of two in a row in last 5bases of handle
            if re.search('('+repeatUnit+'){2,'+str(len(self.sequence))+'}',self.sequence[-5:]):
               #return repeatUnit+' repeated twice in last five'
               self.logstring += '| homopolymer fail 2x'
               return '2x homopolymer in -5'
            # look for repeated bases in stretches of >=three in a row in full sequence 
            if re.search('('+repeatUnit+'){3,'+str(len(self.sequence))+'}',self.sequence):
               #return repeatUnit+' repeated three times'
               self.logstring += '| homopolymer fail 3x'
               return '3x homopolymer'
    
        # same for 2 base combinations
        for repeatUnit in twoBaseCombos:
           if re.search('('+repeatUnit+'){2,'+str(len(self.sequence))+'}',self.sequence[-5:]):
              #return repeatUnit+' repeated twice in last five'
              self.logstring += '| 2bRepcheck fail 2x'
              return '2x 2brepeat in -5'
           if re.search('('+repeatUnit+'){3,'+str(len(self.sequence))+'}',self.sequence):
              #return repeatUnit+' repeated three times'
              self.logstring += '| 2bRepcheck fail 3x'
              return '3x 2brepeat'
    
        # and three base combinations
        for repeatUnit in threeBaseCombos:
           if re.search('('+repeatUnit+'){2,'+str(len(self.sequence))+'}',self.sequence[-5:]):
              #return repeatUnit+' repeated twice in last five'
              self.logstring += '| 3bRepcheck fail 2x'
              return '2x 3brepeat in -5'
           if re.search('('+repeatUnit+'){3,'+str(len(self.sequence))+'}',self.sequence):
              #return repeatUnit+' repeated three times'
              self.logstring += '| 3bRepcheck fail 3x'
              return '3x 3brepeat'
         
        self.logstring += '| Repcheck Pass'
        return False

    def getBlastHits(self):
        """ Function for blasting the handle sequence against the NCBI nt database to identify homologies
        """
        from Bio.Blast import NCBIWWW
        import sys
        import subprocess as sp
        sys.stdout = Unbuffered(sys.stdout)
 
        
        local=True
        if local:
            #localdb='/sw/data/uppnex/blast_databases/nt'
            localdb='/Users/erikborgstrom/localBioInfo/BLASTnt/nt'
            from Bio.Blast.Applications import NcbiblastnCommandline
            from Bio.Blast import NCBIXML
            from cStringIO import StringIO
            import time
            import os
            
            #setting up blast
            database=localdb
            blastsetting = 'strict'
            infile = open('tmp.fa','w')
            infile.write('>tmp\n'+self.sequence+'\n')
            infile.close()
            if blastsetting == 'strict':  cline = NcbiblastnCommandline(query=infile.name, db=database ,evalue=0.001, outfmt=5)#, out='tmp.blastout')
            elif blastsetting == 'sloppy':cline = NcbiblastnCommandline(query=infile.name, db=database ,evalue=0.001, outfmt=5, dust='no',perc_identity=80, task='blastn')#,out='tmp.blastout')
            cline =                               NcbiblastnCommandline(cmd='blastn', outfmt=5, query=infile.name, db=database, gapopen=5, gapextend=2, culling_limit=2)#,out='tmp.blastout')
            print str(cline)
            
            blast_handle = cline.__call__()
            #blastn = sp.Popen(cline.__str__().split(), stdout=sp.PIPE, stderr=sp.PIPE)
            #blastn.wait()
            #stdout, stderr = blastn.communicate()
            #print blastn.returncode
            #print cline.__str__().split()
            #blast_handle = stdout, stderr
            
            #print blast_handle
    
            blast_handle = StringIO(blast_handle[0])
            blast_handle.seek(0)
            #os.remove(infile.name)
        else:
            sys.stdout.write('getting blast hits for handle#'+str(self.id)+'\n')
            result_handle = NCBIWWW.qblast("blastn", "nr", '>tmp\n'+self.sequence,format_type='XML')
            sys.stdout.write('start parsing blast for handle#'+str(self.id)+'\n')
            from cStringIO import StringIO
            blast_handle = StringIO(result_handle.read())
            blast_handle.seek(0)
 
        from Bio.Blast import NCBIXML
        records = NCBIXML.parse(blast_handle)
        hits=0
        for blast_record in records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    perc_identity = float(hsp.identities) 	/	float(hsp.align_length)	*100
                    perc_coverage = float(hsp.align_length)	/	float(blast_record.query_letters)*100
                    if perc_identity >= 90 and perc_coverage >= 90: hits +=1
        self.blastHits = hits

    def checkMatchedByPrimer(self,activePrimers):
        
        from seqdata import revcomp
        from misc import hamming_distance
        
        # go through all active primers in collection
        for seq, name in activePrimers:
           #print ', '.join([name for seq, name in activePrimers])
           
            # go through all the kmers in the sequence of the handle
            for i in range(len(self.sequence)-self.kmerLength):
                primer3prime = seq[-self.kmerLength:]
                prekmer = self.sequence[i:i+self.kmerLength]
      
                # check if the primer three prime end matches the kmer or the revcomp kmer
                for kmer in [prekmer,revcomp(prekmer)]:
                   if primer3prime == kmer: return name
                   distance = hamming_distance( primer3prime, kmer )
                   if distance < self.minHD: return name
                   if distance < self.minHD+1:
                      distance = hamming_distance( primer3prime[-5:], kmer[-5:] )
                      if distance < 1: return name
        return False

    def check3primEnd(self, kmers,revcompkmers,handles,fiveprime=False):
        """ Look at matches between the handle sequence and the kmer collections
        """
        
        from seqdata import revcomp
        from misc import hamming_distance

        # create dictionary of the kmers in the handle sequence
        revcompself = {}
        for i in range(len(self.sequence)-self.kmerLength+1):
           if not fiveprime: kmer = revcomp(self.sequence[i:i+self.kmerLength])
           else:             kmer =         self.sequence[i:i+self.kmerLength]
           try:            revcompself[kmer].append('revcomp self')
           except KeyError:revcompself[kmer] = ['revcomp self']
        
        # check if the active end matches some other kmers in collection
        toCLose = False # initial vaule
        
        if fiveprime: # set the sequence we are checking and if we are looking at the three primer or five prime end of the sequence
           ENDSEQ = revcomp(self.sequence[:self.kmerLength])
           endName = 'first '
        else:
           ENDSEQ = self.sequence[-self.kmerLength:]
           endName = 'last '
        assert len(ENDSEQ) == self.kmerLength, 'Error: the script is trying to check wrong number of end bases\n'
        self.output = '\ngenereated handle#'+str(self.id)+' '+'check '+endName+str(self.kmerLength)+'='+ENDSEQ+':=> ' # give some info for the output
        
        #
        # check if 3'/5' of the handle match any sequence in kmers or revcomp-kmers
        #
        
        # Look for perfect matches of the end sequence to kmer dictionaries
        for dictionary, name in [(revcompself,'self-rc '),(kmers,''),(revcompkmers,'rc ')]:
           if ENDSEQ in dictionary:
              self.output+= endName+str(self.kmerLength)+' ('+ENDSEQ+') perfect '+name+'match to '+' '+str(dictionary[ENDSEQ]);
              self.resonFordeath = name+'kmer match'
              return

        # Look for matches with missmatch
        for dictionary, name in [(revcompself,'self-rc'),(kmers,''),(revcompkmers,'rc '),(revcompself,'self-rc')]:
            for kmer,hits in dictionary.iteritems():
                assert len(kmer) == self.kmerLength, 'Error: kmer of wrong length: '+kmer+' in '+', '.join(hits)
                if kmer.count('N'): continue
                
                # check for distance of full kmer to kmer dictionaries
                distFull = hamming_distance(ENDSEQ,kmer)
                if distFull < self.minHD:
                    toCLose = True;
                    self.output+= str(distFull)+' mm to '+str(hits)+' ('+kmer+') too close,'
                    self.resonFordeath = name+'kmer match'
                    break
                  
                if distFull < self.minHD+1: # if almost to close check last five bases so that we have at least 2 nonmatching bases in this part
                    distLastFive = hamming_distance(ENDSEQ[-5:],kmer[-5:])
                    if distLastFive < 3:
                        toCLose = True;
                        self.output+= str(distLastFive)+' mm in last5 to '+str(hits)+' ('+kmer+') too close,'
                        self.resonFordeath = name+'kmer match  in last5 '
                        break
                
                #if distFull < self.minHD+1 and ENDSEQ[-3] == kmer[-3]:  # looks for uniq three mers skip sthis mostly there are non 4**3 is to few
                #      toCLose = True;
                #      self.output+= ' lastbase(s) identical, '+name+' '+kmer
                #      self.resonFordeath = name+'lastbase(s) identical '
                #      break
                #else:output+= str(dist)+' mm to '+str(hits)+' ok '
            if toCLose: return
    
        #
        # check if 3' bases in handle matches any other 3' in other handles
        #
        for handle2 in handles:
           dist = hamming_distance(ENDSEQ,handle2.sequence[-self.kmerLength:])
           if dist < self.minHD:
              toCLose = True;
              self.output+= str(dist)+' mm to '+str(handle2.id)+'('+handle2.sequence[-self.kmerLength:]+')'+' too close,'
              self.resonFordeath = '3 prime ends match'
              break
           else: self.output+= str(dist)+' mm to '+str(handle2.id)+' ('+handle2.sequence[-self.kmerLength:]+') '+'ok |'+' '
        if toCLose: return   

    def generateAndCheck(self,activePrimers,kmers,revcompkmers,handles):
        """ check if the handle sequence pass filters
        """
        
        #look for repeat structures
        if not self.isRepeat: pass
        else:
          self.dead = True
          return None, self.isRepeat
     
        # check that tm matches set range
        if  self.betterTm >= self.min_tm and self.betterTm <= self.max_tm: self.logstring += '| Tm Pass'
        else:
          self.dead = True
          self.logstring += '| Tm Fail'
          return None, 'tm fail'
     
        #check GC content
        if  self.GCcontent >= self.lower_gc and self.GCcontent <= self.upper_gc: self.logstring += '| GC Pass'
        else:
          self.dead = True
          self.logstring += '| GC Fail'
          return None, 'gc fail'
     
        # look if any of the primers in the database matches towards the handle
        # kmerLength match length and
        matched = self.checkMatchedByPrimer(activePrimers)
        if not matched: pass
        else:
           self.resonFordeath = 'matched by primer ' +matched
        self.check3primEnd(kmers,revcompkmers,handles)

        if self.resonFordeath:
           print str(self.output) + ' ' + str(self.resonFordeath)
           self = None
        else:
           print str(self.output) + ' ' + str(self.resonFordeath) + ' Nice handle!'
           #return self,'good handle'
        return self, 'all ok'

def blastHandle(handle,timeforlastblast):
    #
    # imports
    #
    import time
    import sys
    sys.stdout = Unbuffered(sys.stdout)
    waitTime = 0 
 
    #
    # Balst the handle sequence
    #
    if timeforlastblast and time.time()-timeforlastblast < waitTime:
       sys.stdout.write('blast was too soon waiting ')
       for i in range(waitTime-(time.time()-timeforlastblast)):
          time.sleep(1)
          sys.stdout.write('.')
       sys.stdout.write('\n')
    handle.getBlastHits()
    timeforlastblast = time.time()
    if handle.blastHits:
       print 'handle has more than one good blast hit to nr';
       return None
    return handle

def addHandleToCollection(handle,handles,kmers,revcompkmers,activePrimers,kmerLength):
    from seqdata import revcomp, comp
    handles.append(handle)
    activePrimers.append((handle.sequence,handle.id))
    activePrimers.append((revcomp(handle.sequence),str(handle.id)+'_rc'))
    sequence = handle.sequence
    name = handle.id
    kmers,revcompkmers = toKmers(kmers,revcompkmers,sequence,name,kmerLength)
    return handle,handles,kmers,revcompkmers,activePrimers

def toKmers(kmers,revcompkmers,sequence,name,kmerLength):
    
    from seqdata import revcomp, comp
    
    for i in range(len(sequence)-kmerLength+1):
        kmer = sequence[i:i+kmerLength]
        print kmer
        try: kmers[kmer].append(name)
        except KeyError:kmers[kmer]=[name]
        try: revcompkmers[revcomp(kmer)].append(name)
        except KeyError:revcompkmers[revcomp(kmer)]=[name]
    return kmers,revcompkmers

class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)




