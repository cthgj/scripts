#! /usr/bin/env python

###
###  
###
from _vars import MO_SUBDIR,MYSQL_HOST,MYSQL_DB,MYSQL_USER,MYSQL_PASSWD,GREEN,RED,ensembl_code2name,TREES_REPOSITORIES,PUBLIC_PHYLOMES,TMP,EXTERNAL_ID_DBS
import MySQLdb
import os

class metaPhOrs(object):
  """ Connect to metaPhors database. Return cursor object. """
  def __init__( self, MYSQL_DB, MYSQL_USER=MYSQL_USER, MYSQL_PASSWD=MYSQL_PASSWD, MYSQL_HOST=MYSQL_HOST ):
    self.MYSQL_DB = MYSQL_DB
    self.MYSQL_USER = MYSQL_USER
    self.MYSQL_PASSWD = MYSQL_PASSWD
    self.MYSQL_HOST = MYSQL_HOST
    self.cursor = self.getCursor()
    self.SPECIES,self.code2taxid=self.getSpeciesInfo()
    
  def getCursor( self ):
    cnx=MySQLdb.connect( user=self.MYSQL_USER, passwd=self.MYSQL_PASSWD, host=self.MYSQL_HOST, db=self.MYSQL_DB ) 
    return cnx.cursor()

  def getSpeciesInfo( self ):
    """Return dictionary of all species codes available in meta_orthology as keys 
    and again species code, and their latin species name, and NCBI taxid code as a tuple."""
    tables=[]
    if self.cursor.execute("SHOW tables"):
      for tname, in self.cursor.fetchall():
        try: taxid=int(tname)
        except:continue
        tables.append(taxid)
    sp_names,code2taxid = {},{}
    if self.cursor.execute("SELECT taxid, code, name FROM species" ):
      for taxid, code, name in self.cursor.fetchall():
        if taxid not in tables: 
          if taxid==272634: ###that's really bad, but might work
            sp_names[ taxid ] = ( code, name, 2104 )
            if not code in code2taxid: code2taxid[code]=2104
          continue #if code not in tables: continue 
        sp_names[ taxid ] = ( code, name, taxid )
        if not code in code2taxid: code2taxid[code]=taxid
    return sp_names,code2taxid
  
  def get_external( self,id,single=False ):#,extDB=None,queryDB='' ):
    """Figure external ID and dbs from query id. 
    extDB has to be tuple of db names"""
    if not id.startswith( 'Phy' ): id='Phy'+id
    extIDsData=[ (id,'meta') ]
    cmd0="""SELECT external_id,db FROM ext2protid
            WHERE db IN %s AND protid='%s'
            """ % ( EXTERNAL_ID_DBS,id[3:] )
    if self.cursor.execute( cmd0 ):
      extIDsData=[]
      for tup in self.cursor.fetchall(): extIDsData.append( tup )
    
    if single: return extIDsData[0]
    else: return extIDsData
  def get_external1( self,id,single=False,extDB=('sprot','trembl') ):#,extDB=None,queryDB='' ):
    """Figure external ID and dbs from query id. 
    extDB has to be tuple of db names"""
    if not id.startswith( 'Phy' ): id='Phy'+id
    extIDsData=[ (id,'meta') ]
    cmd0="""SELECT external_id,db FROM ext2protid
            WHERE db IN %s AND protid='%s'
            """ % ( EXTERNAL_ID_DBS,id[3:] )
    if self.cursor.execute( cmd0 ):
      extIDsData=[]
      for tup in self.cursor.fetchall(): extIDsData.append( tup )
    
    if single: return extIDsData[0]
    else: return extIDsData

  def getSeq( self,id,returnSpCode=False ):
    """Return aa seq and taxid of metaids (also PhyID). 
    Return species code instead of taxid if returnSpCode is True"""
    seq=taxid=None
    id=id.split('_')[0]
    if returnSpCode:  cmd="""SELECT seq,code FROM proteins,species WHERE protid='%s' AND proteins.taxid=species.taxid""" % id[3:] 
    else:             cmd="""SELECT seq,taxid FROM proteins WHERE protid='%s'""" % id[3:] 
    if self.cursor.execute(cmd): seq,taxid=self.cursor.fetchone()
    return seq,taxid

  def get_background_color(  self,db_signal,CS_th,orthologs=1,external=True ):
    """Return html code for background color for given db_signal, CS_th and homology_type """
    if db_signal != -1:
      bgcolor = GREEN
      if not external or orthologs: #for phylomeDB cs
        if db_signal<CS_th: bgcolor = RED
        elif not orthologs and db_signal==CS_th and not external: bgcolor = RED #for phylomeDB where CS==CS_th
      else: #for paralogy homology_type = 0
        if db_signal!=0: bgcolor = RED
      return ' bgcolor="%s"' % bgcolor
    return ''

  def get_orthologs_between_2_species( self,taxid1,taxid2,CS_th=0.5,EL_th=1,external_dbs=[],homology_type=1,exclude_db_list=[] ):
    """Return all orthology relationships between given two sp."""
    d = {}
    compare = '>='
    if homology_type == 0: compare = '<' #searching for paralogs
    cmd="""SELECT protid1,protid2,CS,trees,phylomes,signals FROM `%s` WHERE taxid2='%s' AND CS %s %s""" % ( taxid1, taxid2, compare, CS_th )
    if self.cursor.execute( cmd ):
      for protid1,protid2,CS,trees,phylomes,signals in self.cursor.fetchall(): 
        if not self._include( phylomes,signals,exclude_db_list,EL_th,external_dbs,homology_type ): continue  #skip results not passing EL criteria etc
        data_tuple=( protid2,CS,trees,phylomes,signals )
        data_list = list( data_tuple ) #external_id2, CS, otrees, trees,  phylomes, ensembl, eggnog, orthoMCL, cog, fo
        try: d[ protid1 ].append( data_list )
        except: d[ protid1 ] = [ data_list, ]
    return d
    
  def get_orthologs_by_id( self, id, CS_th=0.5, external_dbs=[], homology_type=1, EL_th=1, show_structure='1' ):
    """Return 2 objects:
    1. sorted dictionary (by species) of orthologs for given query id.
    2. list of ( species, external_id, external_db, [internal ids] ) if query id is recognised, 
    otherwise return list containing only one element ( query_id ). """
    orthologs = { }; compare = '>='
    if homology_type == 0: #searching for paralogs
      compare = '<='
      CS_th = 1 - CS_th
    
    query_list = self.get_id( id )
    taxid,query = query_list[:2]

    if taxid: orthologs = self._get_orthologs( taxid,query[3:],compare,CS_th,external_dbs,EL_th,homology_type,show_structure )

    return orthologs, query_list

  def  _get_orthologs(  self,taxid,protid,compare,CS_th,external_dbs,EL_th,homology_type,show_structure ):
    """Internal function, execute db homology_type and return orthology dict"""
    exclude_db_list = [ ]
    orthologs = { }; coorthologs = { }
    id2_list = [ ]
    if taxid in self.SPECIES.keys():
      cmd = """SELECT taxid2,protid2,CS,trees,phylomes,signals FROM `%s` WHERE protid1='%s' AND CS %s %s""" % ( taxid,protid,compare,CS_th )                                   
      if self.cursor.execute( cmd ):
        #SAVE HOMOLOGY DATA INTO DICTIONARY
        for taxid2,protid2,CS,trees,phylomes,signals in self.cursor.fetchall(): # taxid2,external_id1,external_id2,CS,trees,phylomes,extSignals 
          if not self._include( phylomes,signals,exclude_db_list,EL_th,external_dbs,homology_type ): continue #skip results not passing EL criteria etc
          data_list = [ protid2,CS,trees,phylomes,signals, None ]
          id2_list.append( protid2 )
          try: orthologs[ taxid2 ].append( data_list )
          except: orthologs[ taxid2 ] = [ data_list ]
        #ADD CO-HOMOLOGY DATA IF REQUESTED
        if show_structure == '1' and len( id2_list ) > 1:
          cmd = """SELECT protid1, taxid2, protid2, CS FROM `%s` 
                        WHERE protid2 IN %s AND CS %s %s
                        """ % ( taxid, tuple( id2_list ), compare, CS_th )
          if self.cursor.execute( cmd ):
            for protid1, taxid2, protid2, CS in self.cursor.fetchall():
              if protid1 == protid: continue
              co_orthology_tuple = ( protid1, CS ) 
              if coorthologs.has_key( taxid2 ):
                try: coorthologs[ taxid2 ][ protid2 ].append( co_orthology_tuple )
                except: coorthologs[ taxid2 ][ protid2 ] = [  co_orthology_tuple ]
              else: coorthologs[ taxid2 ] = { protid2: [ co_orthology_tuple ] }
        #COMBINE HOMOLOGY AND CO-HOMOLOGY DICTIONARIES    
        for taxid2 in orthologs:
          i = -1
          for datalist in orthologs[ taxid2 ]:
            i += 1
            protid2 = datalist[0]
            try: data = coorthologs[ taxid2 ][ protid2 ]; data.sort(); datalist[-1] = data
            except: continue 
            orthologs[ taxid2 ][ i ] = datalist

    return orthologs
      
  def _include( self,phylomes,signals,exclude_db_list,EL_th,external_dbs,homology_type ):
    """Return true if signals passed by external_data_list are following criteria in exclude_db_list and evidence_level"""
    _include = False
    el = i = 0
    external_data_list=[ phylomes ]# phylomes
    for signal in signals:
      try: external_data_list.append( int(signal) )
      except: external_data_list.append( -1 )
    #CHECK IF GIVEN DATABASE CONFIRM PREDICTION
    for db in TREES_REPOSITORIES:
      if db in external_dbs:
        if homology_type == 1:  #for otrthology predictions
          if external_data_list[ i ] <1: return False # reject prediction if it's not confirmed by particular db, or phylomes no < 1
        else: #for paralogy predictions
          if db=='phylomes': 
            if not external_data_list[ i ]: return False
          elif external_data_list[ i ] != homology_type: return False
      i += 1
      
    # CHECK EVIDENCE LEVEL AND exclude_db_list criterium
    for db,signal in zip( TREES_REPOSITORIES, external_data_list ): #TREES_REPOSITORIES keep db names, external_data_list keep signal from each db
      if db not in exclude_db_list:
        if db=='phylomes': #phylomes can be 0, 1, and higher number
          if signal >= 1: #HOW TO INCLUDE PHYLOMEDB CS HERE??
            _include = True
            el += signal
            
        else: #signals from external dbs can be -1-no info, 0-paralogs, 1-orthologs
          if signal == homology_type:
            _include = True
            el += 1
    
    #DON'T INCLUDE IF EL_TH > that no of dbs giving positive signals
    if el < EL_th: _include = False
        
    return _include

  def get_multi_orthologs( self, queries, species, CS_th=0.5, EL_th=1, external_dbs=[], homology_type=1, exclude_db_list=[] ):
    """Return dictionary of orthology predictions for multiple queries"""
    compare = '>='; queries_str = ''; species_str = ''
    taxids=[]
    if species:
      for taxid in species:
        if type(taxid)!=int: 
          try:taxid=self.code2taxid[taxid]
          except:continue
        taxids.append(taxid)
    if homology_type == 0: compare = '<' #searching for paralogs
    taxid = None; sp_err = 0; queries_dict = {}; orthologs_multi = {}; species_with_homologs=[ ]
    ###get query species from queries list
    for q in queries:
      query_list = self.get_id( q )
      queries_dict[ q ] = query_list
      queries_str += " '%s'," % query_list[1][3:]
      if not taxid: taxid = query_list[0]
      else:
        new_taxid = query_list[0]
        if new_taxid != taxid: sp_err += 1
        
    if not taxid: return orthologs_multi, sp_err, queries_dict, species_with_homologs # return non if species not recognized
    #if not species: species=None #avoid filtering of taxid2 when all species selected
    queries_str = queries_str[:-1]
    cmd = """SELECT taxid2,protid1,protid2,CS,phylomes,signals FROM `%s` WHERE protid1 IN (%s) AND CS%s%s""" % ( taxid, queries_str, compare, CS_th )
    if self.cursor.execute( cmd ): # taxid2,external_id1,external_id2, CS, trees, signals
      for taxid2,protid1,protid2,CS,phylomes,signals in self.cursor.fetchall():
        if taxids and taxid and taxid2 not in taxids: continue
        if not self._include( phylomes,signals,exclude_db_list,EL_th,external_dbs,homology_type ): continue  #return False if prediction didn't pass
        protid1,protid2='Phy'+protid1,'Phy'+protid2
        data = ( protid2,CS )
        if orthologs_multi.has_key( protid1 ):
          try:  orthologs_multi[protid1][taxid2].append( data )
          except: orthologs_multi[protid1][taxid2] = [ data, ]
        else:  orthologs_multi[protid1] = { taxid2 : [ data, ] }
        if taxid2 not in species_with_homologs: species_with_homologs.append( taxid2 )  
    return orthologs_multi,sp_err,queries_dict,species_with_homologs
    
  def get_id( self,id,extDB=True ):
    """Check if input is internal id, if not get internal id and return a list of those.
    If neither phylomeDB, nor conversion is found return monit about wrong ID.
    Return list of ( species, external_id, external_db ) if query id is recognised, 
    otherwise return list containing only one element ( query_id ). """
    for el in [ '\t', '\r', ' ', '\s', '\n' ]: id = id.replace( el, '' )
    #for phylomeDB_uniq ids
    self.cursor = self.getCursor()
    external_ids = [ ]
    if id.startswith( ('Phy','Mpo') ): protid,taxid=self._getInternalID( id )
    else: protid,taxid=self._getByExtID( id ) #for external ids

    if extDB: external_ids=self.get_external( protid,False )
    
    if taxid in self.SPECIES.keys(): return taxid,protid,external_ids
    else: return None,protid,external_ids #or no sp if wrong species

  def _getInternalID( self,id ):
    """Get PhyID or metaID with sp code for given query"""
    taxid=None
    protid=id.split('_')[0]
    if self.cursor.execute( """SELECT taxid FROM proteins WHERE protid='%s'""" % protid[3:] ): 
      taxid= self.cursor.fetchone()[0]
      if taxid==272634: taxid=2104
    
    return protid,taxid

  def _getByExtID( self,id ):
    """Get internalID and sp for externalID"""
    taxid=None
    cmd1="""SELECT DISTINCT ext2protid.protid,taxid FROM ext2protid,proteins WHERE db IN ( 'sprot','trembl','Ensembl' ) AND external_id='%s' 
            AND ext2protid.protid=proteins.protid""" % ( id, )
    cmd3="""SELECT DISTINCT ext2protid.protid,taxid FROM ext2protid,proteins WHERE external_id='%s' AND ext2protid.protid=proteins.protid""" % ( id, )
    cmd2="""SELECT id_conversion.protid,taxid FROM id_conversion,proteins WHERE external_id='%s' AND id_conversion.protid=proteins.protid""" % ( id, )
    if self.cursor.execute( cmd1 ):   id,taxid=self.cursor.fetchone()
    elif self.cursor.execute( cmd2 ): id,taxid=self.cursor.fetchone()
    elif self.cursor.execute( cmd3 ): id,taxid=self.cursor.fetchone()
    if taxid: id='Phy'+id
    return id,taxid
  
