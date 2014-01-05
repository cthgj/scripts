### SERVER SPECIFIC PART!!! PAY ATTENTION TO ALTER THIS
HOME = '/home/lpryszcz/'
HOST = '/home/services/web/orthology.phylomedb.org/'
BLAST_DB_PATH = HOME + 'db_201009/'
BLAST_BIN_PATH = HOST + 'wsgi/blastall -p blastp'
BLAST_DEFAULT_FLAGS = "-FT -m9 -e 1e-25 -b 500 -t0 -a 2"  # -a parameter set number of cores for blast search, maybe 3?
TMP = HOST + 'tmp/'
MO_SUBDIR = ''
MYSQL_DB = 'metaPhOrs_201011'
MYSQL_HOST = 'localhost'#'neptune'
MYSQL_USER = 'script'
MYSQL_PASSWD = 'python'
main_http='http://orthology.phylomedb.org/'

# HTML COLORS CODES
WHITE1 = '#F8F8F8' #almost white #white
GREY1 = '#E0E0E0'  #light grey
GREEN = '#d4ff81' # '#66FF00' 
RED = '#ff3c3c'# '#FF0000' 
TABLE_HEADER_COLOR = '#E0E0E0'  #light grey

TREES_REPOSITORIES=( 'phylomes', 'ensembl', 'eggnog',  'orthoMCL', 'cog', 'fo', 'treefam' )
TREES_REPOSITORIES_NAMES=( 'PhylomeDB', 'Ens', 'Egg', 'Ort', 'COG', 'FO', 'TF' )
TREES_REPOSITORIES_FULLNAMES= ('PhylomeDB', 'Ensembl', 'EggNOG', 'OrthoMCL', 'COG', 'Fungal Orthogroups', 'TreeFAM' )
EXTERNAL_TREES_REPOSITORIES =TREES_REPOSITORIES[1:]
EXTERNAL_REPOSITORIES_NAMES=TREES_REPOSITORIES_NAMES[1:]
EXTERNAL_ID_DBS=( 'sprot','trembl' ) # 'ensembl',

PUBLIC_PHYLOMES=[ 1,3,4,5,7,16,18,19,20,21,23,26,28 ]

DEAFULT_SESSION = {
  'query': 'Phy0007XA1_HUMAN', #'P47058',
  'homology_type': 1,
  'ext_signals': None, 
  'external_dbs': '', 
  'CS_th': 0.5,
  'EL_th': 1,
  'show_structure': '1',
  'format': 'html', 
  'external_dbs': [],
  'coortologs_limit': 1,
  'spName_maxLen': 30,
  }

CONV_DICT = {
         -1: '-',
         0: 'paralog', 
         1: 'ortholog', 
         }

db_link = { 
        'uniprot' : 'http://www.uniprot.org/uniprot/', 
        'sprot':'http://www.uniprot.org/uniprot/', 
        'trembl':'http://www.uniprot.org/uniprot/', 
        'ensembl': 'http://www.ensembl.org/Homo_sapiens/Search/Results?species=all;idx=;q=', 
        'Ensembl': 'http://www.ensembl.org/Homo_sapiens/Search/Results?species=all;idx=;q=', 
        'phylome': 'http://phylomedb.org/?q=node/3&seqid=', 
        'meta': 'http://orthology.phylomedb.org/?q=seqInfo&seqid=', 
	}

ensembl_code2name = {
  'CPO': 'Cavia_porcellus', 
  'OAN': 'Ornithorhynchus_anatinus', 
  'MIC': 'Microcebus_murinus', 
  'Rno': 'Rattus_norvegicus', 
  'VPA': 'Vicugna_pacos', 
  'CHO': 'Choloepus_hoffmanni', 
  'GAC': 'Gasterosteus_aculeatus', 
  'LAF': 'Loxodonta_africana', 
  'EEU': 'Erinaceus_europaeus', 
  'TSY': 'Tarsius_syrichta', 
  'PCA': 'Procavia_capensis', 
  'OCU': 'Oryctolagus_cuniculus', 
  'MEU': 'Macropus_eugenii', 
  'TTR': 'Tursiops_truncatus', 
  'PTR': 'Pan_troglodytes', 
  'CJA': 'Callithrix_jacchus', 
  'RNO': 'Rattus_norvegicus', 
  'BTA': 'Bos_taurus', 
  'TBE': 'Tupaia_belangeri', 
  'SSC': 'Sus_scrofa', 
  'MOD': 'Monodelphis_domestica', 
  'ACA': 'Anolis_carolinensis', 
  'MUS': 'Mus_musculus', 
  'ETE': 'Echinops_telfairi', 
  'DOR': 'Dipodomys_ordii', 
  'ECA': 'Equus_caballus', 
  'TRU': 'Takifugu_rubripes', 
  'PVA': 'Pteropus_vampyrus', 
  'GGO': 'Gorilla_gorilla', 
  'DNO': 'Dasypus_novemcinctus', 
  'SAR': 'Sorex_araneus', 
  'MMU': 'Macaca_mulatta', 
  'XET': 'Xenopus_tropicalis', 
  'ORL': 'Oryzias_latipes', 
  'CIN': 'Ciona_intestinalis', 
  'CAF': 'Canis_familiaris', 
  'TGU': 'Taeniopygia_guttata', 
  'GAL': 'Gallus_gallus', 
  'DAR': 'Danio_rerio', 
  'PPY': 'Pongo_pygmaeus', 
  'STO': 'Spermophilus_tridecemlineatus', 
  'MLU': 'Myotis_lucifugus',
  'TNI': 'Tetraodon_nigroviridis', 
  'OPR': 'Ochotona_princeps', 
  'OGA': 'Otolemur_garnettii'
  }
