#!/usr/bin/python

#initialize the metaPhOrs connection
import dbClient
import sys, getopt, string, re

## Get options
try:
    options, xarguments = getopt.getopt(sys.argv[1:],
    'ha', ['species=','file=', '--view', 'version'])
except getopt.error:
    print 'Options Error'
    sys.exit(0)
for a in options[:]:
    if a[0] == '--file' and a[1] != '':
        filename= a[1]
        options.remove(a)
        
    elif a[0] == '--file' and a[1] == '':
        print '--file expects an argument'
        sys.exit(0)
    if a[0] == '--species' and a[1] != '':
        spec= a[1]
        options.remove(a)
        
    elif a[0] == '--species' and a[1] == '':
        print '--species expects an argument'
        sys.exit(0)
    if a[0] == '--all' :
        get_all=1

m=dbClient.metaPhOrs('metaPhOrs_201011','metapi','python','cgenomics.crg.es')

## Read names 
names={}
f=open(filename,'r')
for line in f:
    names[line.strip()]=1;
if spec=='human':
    species=('DROME','CAEEL','MOUSE','YEAST')
    taxids={7227:'fly',6239:'worm',10090:'mouse',4932:'yeast',}

elif spec=='fly':
    species=('HUMAN','CAEEL','MOUSE','YEAST')
    taxids={9606:'human',6239:'worm',10090:'mouse',4932:'yeast',}
elif spec=='mouse':
    species=('DROME','CAEEL','HUMAN','YEAST')
    taxids={7227:'fly',6239:'worm',9606:'human',4932:'yeast',}
elif spec=='worm':
    species=('DROME','HUMAN','MOUSE','YEAST')
    taxids={7227:'fly',9606:'human',10090:'mouse',4932:'yeast',}
elif spec=='yeast':
    species=('DROME','CAEEL','MOUSE','HUMAN')
    taxids={7227:'fly',6239:'worm',10090:'mouse',9606:'human',}
elif spec=='meta':
    species=('MACMU','RAT','BOVIN','PANTR')
    taxids={7227:'fly',6239:'worm',10090:'mouse',9606:'human',}
elif spec=='all':
    species=('*')
else:
    sys.stderr.write("ERROR: --species must be on of human, fly, yeast, mouse, worm, meta, all\n")
    sys.exit(0)


taxids={9606:'human',7227:'fly',6239:'worm',10090:'mouse',4932:'yeast',10116:'rat',9913:'cow',9598:'chimp',9544:'macaque'}
#print "\t%s" % '\t'.join(map(str,species))

meta_names={};
netw_names={};
num=len(names);
k=1
for name in names:
    sys.stderr.write("Getting names : %s of %s\r" % (k,num)),
    k=k+1
    meta_name=m.get_id(name)[1].strip();
    meta_names[name.strip()]=meta_name;
    netw_names[meta_name]=name.strip();

sys.stderr.write("\n")

    
sys.stderr.write("Getting orths...")
print species

orths=m.get_multi_orthologs(netw_names.keys(),species ) #you can specify CS,EL as well and so on
sys.stderr.write(" Done\n")
print orths
#orths is :
#({'Phy0008C3Q': {9544L: [('Phy000ALW1', 1.0)], 9913L: [('Phy00022YY', 1.0), ('Phy001Q8LN', 1.0), ('PhyZ000BGK', 1.0)], 10116L: [('Phy000CM4N', 1.0)], 9598L: [('PhyZ0034HN', 1.0), ('Phy000BOH5', 1.0)]}, 'Phy00086M9': {9544L: [('PhyZ002DRX', 0.667), ('PhyZ0024O2', 1.0)], 9913L: [('Phy0001ZEQ', 0.83)], 10116L: [('Phy000CAXW', 0.981)], 9598L: [('Phy000BINY', 0.833)]}, 'Phy0008D13': {9544L: [('Phy000ACVW', 1.0), ('PhyZ002113', 1.0)], 9913L: [('Phy0001WG2', 0.979), ('Phy001Q8LO', 1.0), ('PhyZ000AJL', 1.0)], 10116L: [('Phy000CESU', 0.981)], 9598L: [('Phy000BPGJ', 1.0)]}}, 0, {'Phy0008C3Q': (9606L, 'Phy0008C3Q', [('P31946', 'sprot')]), 'Phy0008D13': (9606L, 'Phy0008D13', [('Q04917', 'sprot')]), 'Phy00086M9': (9606L, 'Phy00086M9', [('D3DTH5', 'trembl'), ('P62258', 'sprot')])}, [9544L, 9913L, 9598L, 10116L])




ddd=orths[0];
n=1;
#print orths
for name in names:
    sys.stderr.write("%s of %s\r" % (n,num)),
    n=n+1;
#    try:
    if ddd.has_key(meta_names[name.strip()]):
        new_dict=ddd[meta_names[name.strip()]]
        #sys.stderr.write("Done")
        print "Orig name: %s"%name
        for taxid in new_dict.keys():
            print "TAX : %s",taxid
            orth=new_dict[taxid]
            for suborth in orth:
                ex_names=m.get_external1(suborth[0]);
                #print ex_names
                if ex_names[0][1] != 'meta':
                    print "\t%s"%ex_names[0][0],
        print ""
sys.stderr.write("\n")


