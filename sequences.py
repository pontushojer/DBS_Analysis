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
WFA_H1  = 'CAGTTGATCATCAGCAGGTAATCTGG' #E
WFA_DBS = 'BDHVBDHVBDHVBDHVBDHV'
WFA_H2 = 'CTGTCTCTTATACACATCTCATGAGAACGTCGTTGACGATGGACAGTTCCAAGAGGTCATG' #H1691'+H5+TES
WFA_H3 = 'TAGGACCAGCGTCTCAGTAGAGATGTGTATAAGAGACAG' #H43283'G+TES
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

def main():
    
    import seqdata
    
    print '#  '
    print '#  The expected layout of inserts should be:'
    print '#  '
    print '#  H1-DBS-revcomp(H2)-someDNA-revcomp(H3)'
    print '#  '
    print '#  Using the currently defined sequences this should be:'
    print '#  '+H1+'-'+DBS+'-'+seqdata.revcomp(H2)+'-someDNA-'+seqdata.revcomp(H3)
    print '#  '
    print '#  '+'With illumina handles this will be:'
    print '#  '+ILLI5+'-'+H1+'-'+DBS+'-'+seqdata.revcomp(H2)+'-someDNA-'+seqdata.revcomp(H3)+'-'+ILLI7
    print '#  '+'or if ligated the other direction might also occur:'
    print '#  '+ILLI5+'-'+H3+'-someDNA-'+H2+'-'+seqdata.revcomp(DBS)+'-'+seqdata.revcomp(H1)+'-'+ILLI7

if __name__ == "__main__": main()