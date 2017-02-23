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
