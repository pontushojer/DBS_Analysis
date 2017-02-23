class AnalysisFolder(object):
    """This class represent the analysis output folder,
    it hold the structure for it and track the files within
    it is also the container that enables other info such as settings to be sent around among functions in a controlled way"""
    
    def __init__(self, path, logfile=None):
        """ initiates the object and sets the paths to all files relative to commandline input
        also creates the objects database and settings which will be used extensively during the analysis
        """

        #
        # imports
        #
        import os
        import sqlite3
        from dbs_analysis.metadata.database import Database
        from dbs_analysis.metadata.settings import Settings
        from dbs_analysis.metadata.results  import Results

        # Folders
        self.path = path
        self.logpath = self.path+'/logfiles'
        self.rawdataPath = self.path+'/rawData'
        self.dataPath = self.path+'/data'# output/input files
        self.temp = self.path+'/temp'# temp
        
        self.folders = [self.path,self.logpath,self.rawdataPath,self.dataPath,self.temp]

        #Files
        self.databaseFileName  = self.path+'/database.db'
        self.report = self.path+'/report'# report
        self.fastq_outfile1 = self.dataPath+'/inserts.r1.fastq'
        self.fastq_outfile2 = self.dataPath+'/inserts.r2.fastq'
        self.fastq_outfile3 = self.dataPath+'/inserts.singlets.fastq'
        self.coloredReadMasking = self.dataPath+'/coloredReadMasking'
        self.dbsfastq = self.dataPath+'/rawBarcodeSequencesSortedByAbundance.fq'
        
        self.filenames = [self.databaseFileName,self.report]
        
        # flags
        self.database_in_temp = False
        
        # objects
        self.database  = Database(self.databaseFileName, self)
        self.settings = Settings(self)
        self.results = Results(self)
        
        # if database already exists load the info from it
        if os.path.exists(self.databaseFileName):
            try:
                self.settings.loadFromDb()
                self.results.loadFromDb()
            except sqlite3.OperationalError: pass
        
        if type(self.settings.debug) != bool and type(self.settings.debug) != int: self.settings.debug = eval(self.settings.debug)

    def copy_to_temp(self):
        
        import os
        import shutil
        import sys
        
        pid = os.getpid()
        
        if not os.path.exists(self.settings.temp): sys.stderr.write( '# ERROR: the folder '+self.settings.temp+' does not exist...\n' )
        else:
            
            if not self.settings.temp[-1] == '/': self.settings.temp+='/'
            self.database_original_name = self.databaseFileName.replace('//','/')
            sys.stderr.write('copying database ('+self.database_original_name+')  to temp location '+self.settings.temp+'\n')
            shutil.copy2(self.database_original_name,self.settings.temp+'/database.'+str(pid)+'.db'.replace('//','/'))
            sys.stderr.write('done.\n')
            self.databaseFileName = self.settings.temp+'/database.'+str(pid)+'.db'.replace('//','/')
            self.database  = Database(self.databaseFileName, self)
            self.database_in_temp = True

    def copy_from_temp(self):

        import os
        import shutil
        import sys
        
        if self.database_in_temp:
            sys.stderr.write('copying database ('+self.database.path+')  back to original location ('+self.database_original_name+')...\n')
            shutil.copy2(self.database.path,self.database_original_name)
            sys.stderr.write('done.\n')

    def create(self, ):
        """ This functions creates the folder-structure and the database """

        import os
        for folder in self.folders:
            if not os.path.exists(folder):
                os.mkdir(folder)
        
        self.database.create()
        
    def checkIntegrity(self, ):
        """ A function that checks if all folders are in place """
        
        import os
        
        if not os.path.isdir(self.path): return 'FAIL: The folder does not exist.'
        
        for folder in self.folders:
            if not os.path.isdir(folder): return 'FAIL: The folder structure is broken.'
        
        return 'PASS'
    
    def readindexTsv(self, ):
        """ Function that reads a TSV file and returns the entries as list of dictionaries
        each dictrionary in the list will have the following keys corresponding to the columns in the tsv_file:
        indexReference_name
        index_sequences
        """
        # read tsv file save in dictionary
        
        import sys
        import os
        
        if not self.settings.IndexReferenceTsv:
            msg = 'WARNING: no indexreference file is supplied in the settings please run dbs_change_settings and set the appropriate variable.\n'
            msg+= '         (If run is not indexed you can safely ignore this message).\n'
            sys.stderr.write(msg)
            return
        
        if not os.path.exists(self.settings.IndexReferenceTsv):
            msg = 'ERROR: cannot find the specified indexreference file, please run dbs_change_settings and change the setting.\n'
            sys.stderr.write(msg)
            sys.exit()
    
        #
        # get index tsv file name from settings object and open the file
        #
        index_file = open(self.settings.IndexReferenceTsv)
        
        #
        # create the dictionaries where we are going to add the data
        #
        self.individual_id_fasta_sequences_by_id = {}
        self.individual_id_fasta_sequences_by_sequence = {}
        
        #
        # read each line in the file
        #
        for line in index_file:
            
            # split the columns to two variables
            index_name, index_sequences = line.rstrip().split('\t')
            
            # save in dictionaries
            self.individual_id_fasta_sequences_by_id[index_name] = index_sequences 
            self.individual_id_fasta_sequences_by_sequence[index_sequences] = index_name
        
        # end by closing the input file
        index_file.close()
        