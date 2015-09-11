"""sequences.py hold the handle sequences used for the defined layout
change sequences here and they will be used throughut the analysis
"""

H1  = 'GGAGCCATTAAGTCGTAGCT' #H770
DBS = 'BDHVBDHVBDHVBDHVBDHV'
H2 = 'GACAGTTCCAAGAGGTCATG' #H1691
H3 = 'TAGGACCAGCGTCTCAGTAT' #H4328
ILLI5 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLI7 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'

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