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
            if int(binary[-9]):	self.primaryalignment = False;#self.output += 'not primary alignment\t'
            else:			self.primaryalignment = True
            if int(binary[-10]):	self.passfilter=False;#self.output += 'QC failed\t'
            else:			self.passfilter=True
            if int(binary[-11]):	self.pcrduplicate=True;#self.output += 'is PCR duplicate'
            else:			self.pcrduplicate=False;
        except ValueError:
            output += '.'
        return None
