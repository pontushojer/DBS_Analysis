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

def loadBEDfile(filename):
    """ Function that reads a BED file and returns the entries as list of dictionaries
    each dictrionary in the list will have the following keys corresponding to the columns in the bedfile:
        reference_name
        start_position
        end_position
        entry_name
        value
        strand
    """

    bedDictionary = []    
    
    try:
        bed_file = open(filename)
    
        for line in bed_file:
            reference_name, start_position, end_position, entry_name, value, strand = line.rstrip().split('\t')
            bedDictionary.append( {'reference_name':reference_name, 'start_position':int(start_position), 'end_position':int(end_position), 'entry_name':entry_name, 'value':value, 'strand':strand} )
    
    except TypeError as error: # this is if bedfile is None, this should be changed later on so that bedfile list is a clusterwide template copy pasted to each cluster by eg shutil deepcopy
        pass

    return bedDictionary 
   
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

