#
# This file hold some misc functions such as calculating edit distance between strings or reformatting numbers and strings
# see each function for a more detailed description
#

def bufcount(filename):
    """ returns the number of lines in a file
    """
    
    #
    # open the file and check if it's compressed or not
    #
    import gzip
    if filename.split('.')[-1] in ['gz','gzip']: f = gzip.open(filename)
    else: f = open(filename)
    
    #
    # set initial values
    #
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    
    #
    # do the actual counting and return the number of new line characters found
    #
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
        f.close
    return lines

def hamming_distance(s1, s2):
    """ function that compares two strings and returns the number of no matching characters, taken from https://en.wikipedia.org/wiki/Hamming_distance
    """
    assert len(s1) == len(s2), 'Error: '+str(len(s1)) + ' != ' + str(len(s2))+' '+s1+' '+s2
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein(s1, s2):
    """ levenshtein distance, originally taken from somewhere on the web and then edited
    """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if not s1:
        return len(s2)
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            if c1 == 'N' or c2 == 'N': substitutions -= 1 #if N then no mismatch
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def percentage(count,total):
    """ function that calculate percentage from two numbers total and count and returns the result as a float
    if not able to calculate the percentage due to zerodivision or similar the string "NA" is returned
    """
    if type(None) in [type(count),type(total)]: return 'NA'
    if str in [type(count),type(total)]: return 'NA'
    if 'NA' in [total,count]: return 'NA'
    if float(total) <=0.0: return 'NA'
    #return round(float(count) / float(total),4)
    return round(100* float(count) / float(total),2)

def thousandString(string):
    """ takes a number and returns it formatted as a string with a white space every third number
    """
    if type(string) == type(None): return 'NA'
    if type(string) != str: string = str(int(round(float(string),0)))
    outstring = ''
    for i in range(len(string)):
        outstring += string[-(i+1)]
        if (i+1)%3 == 0: outstring += ' '
    return outstring[::-1]

def sorted_nicely( l ): # funtion "stolen from the internet"
    """ Sort the given iterable in the way that humans expect.""" 
    import re
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Progress():
    """ a progress meter
    """
    
    import sys

    def __init__(self,total, verb='full', logfile=sys.stderr, unit='reads' ,mem=False, printint=0):
        import time
        self.total = total
        self.current = 0
        self.type = verb
        self.logfile = logfile
        self.ltime = time.time()
        self.lcurrent = self.current
        self.lpercentage = 0
        if verb == 'full': self.printint = 5
        elif verb == 'minimal':self.printint = 5
        self.unit = unit
        self.mem = mem
        if printint: self.printint = printint

    def __enter__(self):
        if self.type == 'minimal': self.logfile.write('0%                 50%                 100%\n')
        #                                              ....................................................................................

    def update(self):
        import time
        self.current += 1
        self.percentage = int(round(100*float(self.current)/self.total))
        if self.percentage % self.printint == 0 and self.percentage != self.lpercentage:
            self.stf=int(round((self.total-self.current)/((self.current-self.lcurrent)/(time.time()-self.ltime))))
            if self.type == 'full' and self.logfile: self.logfile.write(
                '#Progress => '+str(self.percentage)+'%, '+
                str( round((self.current-self.lcurrent)/(time.time()-self.ltime),2) )+' '+self.unit+'/second, '+
                time.strftime("%A, %d %b %Y %H:%M:%S",time.localtime())+
                ', left: '+str(self.stf/60/60)+'h '+str(self.stf/60%60)+'min '+str(self.stf%60)+'s')
            if self.mem:
                import resource
                total_memory_used = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
                this_process_memory_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                if total_memory_used/1024/1024 > 1024:
                    self.logfile.write(', using '+str(round(float(total_memory_used)/1024/1024/1024,2))+' ('+str(round(float(this_process_memory_used)/1024/1024/1024,2))+') GB.\n')
                elif total_memory_used/1024 > 1024:
                    self.logfile.write(', using '+str(round(float(total_memory_used)/1024/1024,2))+' ('+str(round(float(this_process_memory_used)/1024/1024,2))+') MB.\n')
                else:
                    self.logfile.write(', using '+str(round(float(total_memory_used)/1024,2))+' ('+str(round(float(this_process_memory_used)/1024,2))+') KB.\n')
            else:    self.logfile.write('\n')
            if self.type == 'minimal': self.logfile.write('..')
            self.ltime = time.time()
            self.lcurrent = self.current
            self.lpercentage = self.percentage

    def __exit__(self, *args):
        if self.logfile: self.logfile.write('\n')