class Settings(object,):
    """ The settings objects stores and saves variables in the database in a nice and controlled manner"""
    
    def __init__(self, analysisfolder):
        """ object holding the settings used for each part of the analysis """
        
        #
        # NOTE: All variables needs to be defined in the Defaults,Explenation and Variable sections below!
        #       adding and removing variables can easily be done by modifying these three sections
        #       without need to modify anything else
        #
        
        self.analysisfolder = analysisfolder
        import multiprocessing
        from dbs_analysis import sequences
        
        #
        # the skip list all __dict__ items with the following name will be skipped.
        #
        self.skip_list = ['explenations','defaultValues','isDefault','setTime','analysisfolder','skip_list','TypeOfValues']
        
        #
        # Type of values for all settings
        #
        self.TypeOfValues = {
            'debug':'Bool',
            'uppmaxProject':'string',
            'parallelProcesses':'integer',
            'maxHandleMissMatches':'integer',
            'maxIndividualIndexMissMatches':'integer',
            'barcodeLength':'integer',
            'barcodeMissmatch':'integer',
#            'readsPerClusterCutOff':'integer',
            'bowtie2Reference':'filename',
            'picardPath':'path',
            'mapqCutOff':'integer',
            'minReadDepthForVariantCalling':'integer',
            'minAlleleFreqForHetrozygousVariant':'float',
            'minFreqForSeqPosition':'integer',
            'minPairsPerCluster':'integer',
            'minBasePhredQuality':'integer',
            'targetRegionBed':'filename',
            'IndexReferenceTsv':'filename',
            'temp':'path',
            'port':'string',
            'type':'string',
            'known_hla_types':'filename'
        }
        
        #
        # Default values for all settings
        #
        self.defaultValues = {
            'debug':False,
            'uppmaxProject':'b2014005',
            'parallelProcesses':multiprocessing.cpu_count(),
            'maxHandleMissMatches':0,
            'maxIndividualIndexMissMatches':0,
            'barcodeLength':len(sequences.HLA_DBS),
            #'analysisParts':None,
            'barcodeMissmatch':0,
            #'readsPerUmiCutOff':5,
            #'umiMaxMisMatch':2,
#            'readsPerClusterCutOff':0,
            'bowtie2Reference':None,
            'picardPath':None,
            'mapqCutOff':0,
            'minReadDepthForVariantCalling':10,
            'minAlleleFreqForHetrozygousVariant':0.2,
            'minFreqForSeqPosition':80,
            'minPairsPerCluster':2,
            'minBasePhredQuality':20,
            'targetRegionBed':None,
            'IndexReferenceTsv':None,
            'temp':None,
            'port':'random',
            'type':'HLA',
            'known_hla_types':None
        }
        
        #
        # Explenation strings for all settings
        #
        self.explenations = {
            'debug':'Flag for running the scripts in multiprocessing or as single process run [True/False] (default=False)',
            'uppmaxProject':'Project id used at uppmax for sbatch scripts [bXXXXXXX] (default=b2014005)',
            'parallelProcesses':'Number of process to run when doing multiprocess parts of analysis (defaul=multiprocessing.cpu_count())',
            'barcodeLength':'The length of the bead barcode (default='+str(len(sequences.HLA_DBS))+')',
            #'analysisParts':'Parts of the analysis to run specific for each run.',
            'barcodeMissmatch':'Number of missmatches allowed in the barcode sequence',
            'maxHandleMissMatches':'Number of missmatches allowed in the handle sequence',
            'maxIndividualIndexMissMatches':'Number of missmatches allowed in the individual index sequence',
            #'readsPerUmiCutOff':'Number of reads supporting one UMI for it to passs filters',
            #'umiMaxMisMatch':'Number of missmatches allowed in the UMI sequence',
#            'readsPerClusterCutOff':'Number of reads supporting a barcode sequence cluster for it to passs filters (default = 0)',
            'bowtie2Reference':'path to the bowtie 2 reference index',
            'picardPath':'path to the picard installation to use',
            'mapqCutOff':'filter all reads with mapping quality less than this (default='+str(self.defaultValues['mapqCutOff'])+')',
            'minReadDepthForVariantCalling':'the minimum read depth to call a position as valuable variant calling ( default = 10)',
            'minAlleleFreqForHetrozygousVariant':'the minimum allele frequency of the minor allele at a position to call the position as a hetroxygous variant (0-1, default = 0.2)',
            'minFreqForSeqPosition': 'the minimum freuency of a specific base or insertion at a referrence position to call it as non reference (ie a variant) (0-100 %, default = 80)',
            'minPairsPerCluster':'minimum number of read pairs supporting a cluster for it to be included in analysis (default 2)',
            'minBasePhredQuality':'The minimum base phred quality at the position to call the position as a valuable base (0-40, default = 20)',
            'targetRegionBed':'A bedfile defining regions on the reference that will be used as a targets during analysis such as coverage stats and variant calling (default = None)',
            'IndexReferenceTsv':'A tsv-file of individual index sequnces that will be used as a reference target during analysis of individual index, column 1 should be the name or id and column two should be the sequence (default = None)',
            'temp':'path to temporary folder (used for copying database etc on eg uppmax to increase speed of IO)',
            'port':'web server port for the web interface eg. 5000 (default value = random)',
            'type':'Specify what type of analysis ie what sequences to use to find barcodes etc (default: HLA)',
            'known_hla_types':'fasta file with known HLA types that the consensus sequences will be matched towards (full description in fasta header will be used as identifier for the HLA sequence)'
        }
        
        #
        # containers needed to store information about latest time of seting update to deceide if update is needed
        #
        self.isDefault = {}
        self.setTime = {}
        
        #
        # variables defenition so that they exist in self.__dict__
        #
        self.debug = None
        self.uppmaxProject = None
        self.parallelProcesses = None
        self.maxHandleMissMatches = None
        self.barcodeLength = None
        #self.analysisParts = None
        self.barcodeMissmatch = None
        self.maxIndividualIndexMissMatches = None
        #self.readsPerUmiCutOff = None
        #self.umiMaxMisMatch = None
#        self.readsPerClusterCutOff = None
        self.bowtie2Reference = None
        self.picardPath = None
        self.mapqCutOff = None
        self.minReadDepthForVariantCalling = None
        self.minAlleleFreqForHetrozygousVariant = None
        self.minFreqForSeqPosition = None
        self.minPairsPerCluster = None
        self.minBasePhredQuality = None
        self.targetRegionBed = None
        self.IndexReferenceTsv = None
        self.temp = None
        self.port=None
        self.type=None
        self.known_hla_types=None
        
        # set the default values
        self.setDefaults()

    def setDefaults(self,):
        """ sets the default valeus to the variables using self.__dict__ and matching self.default  """
        for variableName, value in self.defaultValues.iteritems():
            self.__dict__[variableName] = value
            self.isDefault[variableName] = True
            self.setTime[variableName] = None
        return 0

    def loadFromDb(self,):
        """ Load any info from the database settings table and store it this object
        """
      
        #
        # Get the connection
        #
        self.analysisfolder.database.getConnection()

        #
        # Select data
        #
        data = self.analysisfolder.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()

        #
        # Parse data and add to object __dict__
        #
        if data:
            for variableName,default,value,setTime in data:
                self.__dict__[variableName]  = value
                self.isDefault[variableName] = default
                self.setTime[variableName]   = setTime

        #
        # close connection
        #
        self.analysisfolder.database.commitAndClose()

    def setVariable(self,variableName,value):
        """ update the variable value and the set time etc """
        
        import time
        
        # check that the variable name is valid
        assert variableName in self.explenations,'Error: you are trying to set an undefined variable.\n'
        
        # set the values
        self.__dict__[variableName]  = value
        self.isDefault[variableName] = False
        self.setTime[variableName]   = time.time()
        
        return 0

    def saveToDb(self,):
        """ saves the variables in the object to the database table so that they can be loaded in the next executable or process etc...
        """

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
        self.analysisfolder.logfile.write('checking whats in db.\n')
        alreadyInDb = {}
        data = self.analysisfolder.database.c.execute('SELECT variableName,defaultValue,value,setTime FROM settings').fetchall()
        if data:

            for variableName,default,value,setTime in data:
                self.analysisfolder.logfile.write('processing variable '+variableName+'')
                alreadyInDb[variableName] = True

                if variableName in self.__dict__:
                    # check if the variable is the default value or needs update based on time of the set action
                    if default and not self.isDefault[variableName] or setTime < self.setTime[variableName]:
                        
                        # convert any dicts or lists to strings
                        if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
                        
                        # send a message to the log and update database
                        self.analysisfolder.logfile.write(', updating from '+str(value)+' to '+str(self.__dict__[variableName])+', old_setTime '+str(setTime)+' new_setTime '+str(self.setTime[variableName])+'.\n')
                        self.analysisfolder.database.c.execute('UPDATE settings SET defaultValue=?, value=?, setTime=? WHERE variableName=?', (self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName],variableName))
                    
                    else: self.analysisfolder.logfile.write(' no update needed.\n')

        #
        # Add new vars to database
        #
        self.analysisfolder.logfile.write('adding new vars to db:\n')
        for variableName in self.__dict__:
            
            # if not variable but a known container of other types of info such as explenations or settimes
            if variableName in self.skip_list:continue
            
            # new sets ie not updates which were already handled above
            if variableName not in alreadyInDb:
                
                # to str conversion like earlier
                if type(self.__dict__[variableName]) in [dict,list]: self.__dict__[variableName] = str(self.__dict__[variableName])
                
                # update the database
                values = (variableName,self.isDefault[variableName],self.__dict__[variableName],self.setTime[variableName])
                self.analysisfolder.database.c.execute('INSERT INTO settings VALUES (?,?,?,?)', values)
                
                # some info to logfile
                self.analysisfolder.logfile.write('variable '+variableName+' added to db with value '+str(self.__dict__[variableName])+',')
                if self.isDefault[variableName]:self.analysisfolder.logfile.write(' this is the default value.\n')
                else:self.analysisfolder.logfile.write(' non-default value.\n')
            
            else: pass#SEAseqPipeLine.logfile.write('variable\t'+variableName+'\talready in db.\n')

        self.analysisfolder.logfile.write('commiting changes to database.\n')
        self.analysisfolder.database.commitAndClose()
        
        return 0

