#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 8 sept. 2010

@author: benoit
'''


import sys, random, subprocess, Image, os
from optparse import OptionParser


class analyse_clas(object):
    
    def __init__(self, opt, filename):
        
        goldfile = "/home/benoit/Bureau/Guenocheries/ClasChauv/GoldStandard_oldID.txt"
        (mono, multi, distrib, gold, roc, rand) = opt
        
        (list_prot, dict_prot, val) = self.analyse(filename)
        
        if mono:
        	self.writeMono(dict_prot)
        if multi:
        	self.writeMulti(dict_prot)
        if distrib:
        	self.distribClas(val)
        	self.writeDistribClass(dict_prot)
        if gold:
        	self.compare2Gold(dict_prot, list_prot, goldfile)
        if roc:
        	self.writeROCStat(dict_prot, rand)
        
        
        
        
    def writeMulti(self, dict_prot):
        
        list_multi = []
        for prot in dict_prot:
            if dict_prot[prot] > 1:
                list_multi.append(prot)
        
        f = open(filename[:-4]+"multi", "w")
        for prot in list_multi:
            f.write("%s\n"%prot)
        f.close()
        
        
        
    def writeROCStat(self, dict_prot, rand):
    	
    	list_multi = []
        for prot in dict_prot:
            if dict_prot[prot] > 1:
                list_multi.append(prot)
        
        if rand:
        	labels = random.sample(xrange(200), 40)
        else:
        	labels = range(10)
        	labels.extend(range(50,60))
        	labels.extend(range(100,110))
        	labels.extend(range(150,160))
        
        set_labels = set(labels)
        
        tmp = map(lambda x: int(x), list_multi)
        set_multi = set(tmp)
        
        P = 40
       	N = 160
       	T = 200
        TP = len(set_labels.intersection(set_multi))
        FP = len(list_multi) - TP
        FN = P-TP
        TN = T-len(list_multi)-FN
        TPR = float(TP) / P
        FPR = float(FP) / N
        ACC = float((TP + TN)) / (P + N)
        SPC = float(TN) / N
        
        print "file\tTP\tFP\tTN\tFN\tACC\tTPR\tFPR\tSPC"
        print "file\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f" % (TP,FP,TN,FN,ACC,TPR,FPR,SPC)
        
        fn, G = os.path.split(filename)[0]+"/RocStat.txt", os.path.split(filename)[1]
        f=open(fn, "a")
        f.write("%s\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\n" % (G,TP,FP,TN,FN,ACC,TPR,FPR,SPC))
        f.close()
        
        
        
    def writeMono(self, dict_prot):
        
        list_mono = []
        for prot in dict_prot:
            if dict_prot[prot] == 1:
                list_mono.append(prot)
        
        f = open(filename[:-4]+"mono", "w")
        for prot in list_mono:
            f.write("%s\n"%prot)
        f.close()        
        
        
        
    def writeDistribClass(self, dict_prot):
    
        f = open(filename[:-4]+"distribC", "w")
        for prot in dict_prot:
            f.write("%s\t%i\n" % (prot, dict_prot[prot]))
        f.write("\n\n##################\n\n")
        
        val = dict_prot.values()
        list_class = []
        for i in range(max(val)+1):
            list_class.append(val.count(i))
        
        for i in range(len(list_class)):
            f.write("%i\t%i\n" % (i, list_class[i]))
        
        f.close()
        
        
        
    def distribClas(self, val):
        
        f = open("val_tmp.txt", "w")
        for v in val:
            f.write("%.4f\n" % v)
        f.close()
        
        mycmd = """#!/bin/bash
				echo "
				valeurs <- read.table('val_tmp.txt')
				png('hist_scores.png')
				hist(valeurs[['V1']])
				dev.off()
				" | R --slave"""
        try:
            retcode = subprocess.call(mycmd, shell=True)
            #if self.verbose:
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            else:
                print >>sys.stderr, "Child returned", retcode
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e
        
        im = Image.open("hist_scores.png")
        im.show()
        
    
    
    def read_file(self, filename):
    	
        f = file(filename, "r")
        lines = f.readlines()[1:]
        f.close()
        
        list_clas = []
        for i in range(1, len(lines), 2):
            list_clas.append(lines[i].strip().split())
        
        return list_clas
    
    
    
    def readGold(self, arg0):
        
        f = file(arg0, "r")
        lines = f.readlines()
        f.close()
        
        list_gold = map(lambda x: x.strip(), lines[1:])
        
        return list_gold
    
    
    
    def compare2Gold(self, dict_prot, list_prot , goldfile):
        list_gold = self.readGold(goldfile)
        
        list_multi = []
        for prot in dict_prot:
            if dict_prot[prot] > 1:
                list_multi.append(prot)
        
        list_gold2 = []
        list_gold3 = []
        for gold in list_gold:
            if gold in list_multi:
                list_gold2.append(gold)
            if gold in list_prot:
                list_gold3.append(gold)
        
        print "protéines du gold standard dans les classes : %i" % len(list_gold3)
        print "protéines du gold standard multi-classées (%i/%i) :\n%s" % (len(list_gold2), len(list_gold3), ", ".join(list_gold2))
        
    
    
    def analyse(self, filename):
        
        list_clas = self.read_file(filename)
        
        list_prot_all = reduce(lambda x, y: x+y, list_clas)
        list_prot = list(set(list_prot_all))
        #list_prot = list_prot[1:]
        size = len(list_prot_all)
        moy = float(size)/len(list_clas)
        
        dict_prot = {}
        for prot in list_prot:
            dict_prot[prot] = list_prot_all.count(prot)
        
        val = dict_prot.values()
        mono = val.count(1)
        multi = len(dict_prot)-mono
        #maxi = max(val)
        #print maxi
        
        print "nombre de protéines dans les classes : %i" % size
        print "nombre de classes : %i"%len(list_clas)
        print "nombre moyen de protéines par classe : %.3f" % moy
        print "nombre de protéines dans la partition : %i"%len(list_prot)
        print "monoclassées : %i"%mono
        print "multiclassées : %i"%multi
        
        return (list_prot, dict_prot, val)
        
    
    
if __name__ == "__main__":
    
    
    usage = "usage: python %prog [options] filename"
    parser = OptionParser(usage)
    parser.add_option("-m", "--write_mono", action = "store_true", dest="mono", default=False,
                      help="write list of mono-clustered proteins [default: False]")
    parser.add_option("-M", "--write_Multi", action = "store_true", dest="multi", default=False,
                      help="write list of multi-clustered proteins [default: False]")
    parser.add_option("-d", "--distrib", action = "store_true",dest="distrib", default=False,
                      help="display a graph with distribution of classes size and values are written [default: False]")
    parser.add_option("-g", "--gold", action = "store_true", dest="gold", default=False,
                      help="search moonlighting proteins in multi-clutered proteins [default: False]")
    parser.add_option("-r", "--roc", action = "store_true", dest="roc", default=False,
                      help="ONLY FOR SIMULATONS, compute ROC statistiques and write them [default: False]")
    parser.add_option("-R", "--random_nodes", action = "store_true", dest="random", default=False,
                      help="ONLY FOR SIMULATONS, compute ROC statistiques on a random reference and write them [default: False]")
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    
    mono = options.mono
    multi = options.multi
    distrib = options.distrib
    gold = options.gold
    roc = options.roc
    rand = options.random
    
    
    opt = (mono, multi, distrib, gold, roc, rand)
    filename = args[0]
    
    analyse_clas(opt, filename)
    
    
