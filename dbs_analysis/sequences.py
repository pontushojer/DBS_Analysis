"""sequences.py hold the handle sequences used for the defined layout
change sequences here and they will be used throughut the analysis
"""

#################### WFA system ##############################################################################
#H1  = 'CAGTTGATCATCAGCAGGTAATCTGG'#'GGAGCCATTAAGTCGTAGCT' #H770
#DBS = 'BDHVBDHVBDHVBDHVBDHV'
#H2 = 'GACAGTTCCAAGAGGTCATG' #H1691
#H3 = 'TAGGACCAGCGTCTCAGTAT' #H4328
#################### WFA system ##############################################################################

#################### WFA2 system ##############################################################################
import seqdata
WFA_H1  = 'CAGTTGATCATCAGCAGGTAATCTGG' #E
WFA_DBS = 'BDHVBDHVBDHVBDHVBDHV'
#WFA_H2 = 'CTGTCTCTTATACACATCTCATGAGAACGTCGTTGACGATGGACAGTTCCAAGAGGTCATG' #H1691'+H5+TES # Removed for Nextera run
WFA_H2 = 'CTGTCTCTTATACACATCTGACAGTTCCAAGAGGTCATG' # H1691'+TES
WFA_H3 = 'TAGGACCAGCGTCTCAGTAT' # In silico added to 5' end of read 2. Should be isch H4328'
#WFA_H3 = 'TAGGACCAGCGTCTCAGTAGAGATGTGTATAAGAGACAG' #H43283'G+TES # Removed for Nextera run
#################### WFA system ##############################################################################

#################### HLA system ##############################################################################
import seqdata
HLA_H1  = 'ACCGAGTGGTGAGTCATAGT'
HLA_DBS = 'BDVHBDVHBDVHBDVHBDVH'
HLA_H2 = seqdata.revcomp('CTAGCTTCACGAGTTCATCG')
HLA_H3 = 'AGATGGCCGTTATGATAGCG'
#################### HLA system ##############################################################################

#################### Universal ##############################################################################
ILLI5 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLI7 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
IND_HANDLE_1 = 'TTAGTCTCCGACGGCAGGCTTCAAT'
IND_HANDLE_2 = 'ACGCACCCACCGGGACTCAG'
#################### Universal ##############################################################################

def sequence_layout(layout='HLA'):

    if layout == 'HLA':
        H1 = HLA_H1
        H2 = HLA_H2
        H3 = HLA_H3
        DBS= HLA_DBS
    elif layout == 'WFA':
        H1 = WFA_H1
        H2 = WFA_H2
        H3 = WFA_H3
        DBS= WFA_DBS
    else:
        print 'Error: No layout specified.'
        return 1

    import seqdata
    
    output  = '#  \n'
    output += '#  The expected layout of inserts should be:\n'
    output += '#  \n'
    output += '#  H1-DBS-revcomp(H2)-someDNA-revcomp(H3)\n'
    output += '#  \n'
    output += '#  Using the currently defined sequences this should be:\n'
    output += '#  '+H1+'-'+DBS+'-'+seqdata.revcomp(H2)+'-someDNA-'+seqdata.revcomp(H3)+'\n'
    output += '#  '+'\n'
    output += '#  '+'With illumina handles this will be:'+'\n'
    output += '#  '+ILLI5+'-'+H1+'-'+DBS+'-'+seqdata.revcomp(H2)+'-someDNA-'+seqdata.revcomp(H3)+'-'+ILLI7+'\n'
    output += '#  '+'or if ligated the other direction might also occur:'+'\n'
    output += '#  '+ILLI5+'-'+H3+'-someDNA-'+H2+'-'+seqdata.revcomp(DBS)+'-'+seqdata.revcomp(H1)+'-'+ILLI7+'\n'
    
    return output

def main():
    print sequence_layout(layout='HLA')
    print sequence_layout(layout='WFA')


if __name__ == "__main__": main()