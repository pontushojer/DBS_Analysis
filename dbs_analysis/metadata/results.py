class Results(object,):
    
    #
    # this object has the same basic structure as settings see the settings object for comments
    # in the future this might be better handled by inheretance
    #
    
    def __init__(self, analysisfolder):
        """ object holding the results of the analysis """

        self.analysisfolder = analysisfolder
        
        self.defaultValues = {
            'totalReadCount':None,
            'uniqueBarcodeSequences':None,
            'readPairsHasBarcode':None,
            'readPairsAreIlluminaAdapters':None,
            'barcodeClusterCount':None,
            'singeltonBarcodeClusters':None,
            'minR1readLength':None,
            'minR2readLength':None,
            'readsWithDbsPatternMatch':None,
            'constructTypes':None,
            'bt2AlignmentRate':None,
            'alignmentRateQ20':None,
            'alignmentCountQ20':None,
            'alignmentCount':None
        }
        self.explenations = {
            'totalReadCount':'The number of reads totally included in the analysis.',
            'uniqueBarcodeSequences':'The total number of uniqu barcode sequences identified.',
            'readPairsHasBarcode':'Total number of readpairs where a barcode could be identified.',
            'readPairsAreIlluminaAdapters':'Total number of readpairs where a illumina adapter sequence could be identified.',
            'barcodeClusterCount':'Total number of barcode clusters identified.',
            'singeltonBarcodeClusters':'Total number of barcode clusters of only one read pair identified.',
            'minR1readLength':'Minimum read length found in infiles',
            'minR2readLength':'Minimum read length found in infiles',
            'readsWithDbsPatternMatch':'dictrionary holding the counts for the dbs matching statistics of the reapopulation',
            'constructTypes':'dictionary holding infromation of all the different types of constructs found in the read populations',
            'bt2AlignmentRate':'the observed rate of aligned reads',
            'alignmentRateQ20':'rate of SE reads with mappingQ >= 20',
            'alignmentCountQ20':'number of reads with mappingQ >= 20',
            'alignmentCount':'number read pairs in bamfile both mapped and unmapped'
        }
        self.isDefault = {}
        self.setTime = {}
  
        self.totalReadCount = None
        self.uniqueBarcodeSequences = None
        self.readPairsHasBarcode = None
        self.readPairsAreIlluminaAdapters = None
        self.barcodeClusterCount = None
        self.singeltonBarcodeClusters = None
        self.minR1readLength = None
        self.minR2readLength = None
        self.readsWithDbsPatternMatch = None
        self.constructTypes = None
        self.bt2AlignmentRate = None
        self.alignmentRateQ20 = None
        self.alignmentCountQ20 = None
        self.alignmentCount = None
  
        self.setDefaults()

    def setDefaults(self,):
        for resultName, value in self.defaultValues.iteritems():
            self.__dict__[resultName] = value
            self.isDefault[resultName] = True
            self.setTime[resultName] = None
        return 0

    def loadFromDb(self,):
        
        #
        # Get the connection
        #
        self.analysisfolder.database.getConnection()
        
        #
        # Select data
        #
        data = self.analysisfolder.database.c.execute('SELECT resultName,defaultValue,value,setTime FROM results').fetchall()
        
        #
        # Parse data and add to object __dict__
        #
        if data:
            for resultName,default,value,setTime in data:
                self.__dict__[resultName]  = value
                self.isDefault[resultName] = default
                self.setTime[resultName]   = setTime
        
        #
        # close connection
        #
        self.analysisfolder.database.commitAndClose()

    def setResult(self,resultName,value):
        import time
        assert resultName in self.explenations,'Error: you are trying to set a value for an undefined result.\n'
        self.__dict__[resultName]  = value
        self.isDefault[resultName] = False
        self.setTime[resultName]   = time.time()
        return 0

    def saveToDb(self,):
        
        #
        # imports
        #
        import time
        
        #
        # get connection
        #
        self.analysisfolder.database.getConnection()
        
        #
        # Look whats already in database, update it if older or default and set what is not
        #
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('checking results in db.\n')
        alreadyInDb = {}
        data = self.analysisfolder.database.c.execute('SELECT resultName,defaultValue,value,setTime FROM results').fetchall()
        if data:
            for resultName,default,value,setTime in data:
                if self.analysisfolder.logfile: self.analysisfolder.logfile.write('processing result '+resultName+'')
                alreadyInDb[resultName] = True
                
                if resultName in self.__dict__:
                    if default and not self.isDefault[resultName] or setTime < self.setTime[resultName]:
                        if type(self.__dict__[resultName]) in [dict,list]: self.__dict__[resultName] = str(self.__dict__[resultName])
                        if self.analysisfolder.logfile: self.analysisfolder.logfile.write(', updating from '+str(value)+' to '+str(self.__dict__[resultName])+', old_setTime '+str(setTime)+' new_setTime '+str(self.setTime[resultName])+'.\n')
                        self.analysisfolder.database.c.execute('UPDATE results SET defaultValue=?, value=?, setTime=? WHERE resultName=?', (self.isDefault[resultName],self.__dict__[resultName],self.setTime[resultName],resultName))
                    else:
                        if self.analysisfolder.logfile: self.analysisfolder.logfile.write(' no update needed.\n')
        
        #
        # Add new vars to database
        #
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('adding new results to db:\n')
        for resultName in self.__dict__:
            if resultName in ['explenations','defaultValues','isDefault','setTime','analysisfolder']:continue
            if resultName not in alreadyInDb and not self.isDefault[resultName]:
                  if type(self.__dict__[resultName]) in [dict,list]: self.__dict__[resultName] = str(self.__dict__[resultName])
                  values = (resultName,self.isDefault[resultName],self.__dict__[resultName],self.setTime[resultName])
                  self.analysisfolder.database.c.execute('INSERT INTO results VALUES (?,?,?,?)', values)
                  if self.analysisfolder.logfile: self.analysisfolder.logfile.write('result '+resultName+' added to db with value '+str(self.__dict__[resultName])+'\n')
                  #if self.isDefault[resultName]:SEAseqPipeLine.logfile.write(' this is the default value.\n')
                  #else:SEAseqPipeLine.logfile.write(' non-default value.\n')
            else: pass#SEAseqPipeLine.logfile.write('variable\t'+resultName+'\talready in db.\n')
        
        if self.analysisfolder.logfile: self.analysisfolder.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()
        
        return 0

