#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

class CleanGraph(object):
    """
    Nettoie un interactome
    
    A partir d'un graphe, élimine les interactions en double,
    les auto-interactions,
    les conmposantes connexes minoritaires.
    
    Ecrit un nouveau graphe nettoyé
    """
    
    def __init__(self, arg0):
        
        self.cleanGraph(arg0)
        
        
        
    #done
    def readInteractome(self, arg0):
        
        f1 = file(arg0, "r")
        graphe = f1.readlines()
        f1.close()
        graphe = graphe[1:]
        print "nombre d'interactions avant nettoyage : ", len(graphe)
        
        return graphe
    
    
       
    #done
    def computeCC(self, graphe):
        
        graphe = map(lambda x: x.strip(), graphe)
        list_couple = map(lambda x: tuple(x.split('\t')), graphe)
        
        l = []
        for i in range(len(list_couple)):
            l.extend(list(list_couple[i]))
        listeSommets = list(set(l))
        listeSommets2 = listeSommets[:]
        
        dict_compo = {}
        compteur = 0
        for prot in listeSommets:
                compo = set([prot])
                while True:        
                    stable = True
                    for arc in list_couple:
                        if arc[0] in compo and arc[1] not in compo:
                            compo.add(arc[1])
                            del listeSommets[listeSommets.index(arc[1])]
                            stable = False
                        elif arc[1] in compo and arc[0] not in compo: 
                            compo.add(arc[0])
                            del listeSommets[listeSommets.index(arc[0])]
                            stable = False
                    if stable:
                        dict_compo[compteur] = compo
                        compteur = compteur+1
                        break
        
        length = 0
        i = None
        for key in dict_compo.keys():
            if len(dict_compo[key]) > length:
                length = len(dict_compo[key])
                i = key
        
        blackList = list(set(listeSommets2).difference(dict_compo[i]))
        
        print "nombre de protéines à effacer : ", len(blackList)
        #print "liste des protéines à effacer: ", blackList
        
        return blackList
    
    
    #done
    def writeCleanFile(self, arg0, dict_inter):
        
        compt = 0
        for key in dict_inter.keys():
            compt += len(dict_inter[key])
        print "nombre d'interactions après nettoyage : %i"%compt
        
        f = file(arg0[:-3]+"_Clean.gr", "w")
        f.write(str(compt)+"\n")
        sorted_list = sorted(dict_inter.keys())
        for i in sorted_list:
            for elt in dict_inter[i]:
                f.write("%s\t%s\n" %(i, elt))
        f.close()
    
    
    #done
    def cleanGraph(self, arg0):
        
        graphe = self.readInteractome(arg0)
        
        blackList = self.computeCC(graphe)
        i = 0
        state = True
        compt = 0
        while state:
            for j in range(len(blackList)):
                compt = 0
                if blackList[j] in graphe[i]:
                    graphe.remove(graphe[i])
                    compt = 1
                    if i == len(graphe):
                        state = False
                    break
            if compt == 0:
                if i+1 == len(graphe):
                    state = False
                else:
                    i += 1
        
        # élimine les self_loop
        dict_inter = {}
        for i in range(len(graphe)):
            graphe[i]= graphe[i].split('\t')
            graphe[i][1] = graphe[i][1].strip()
            if graphe[i][0] == graphe[i][1]:
                continue
            if graphe[i][0] not in dict_inter.keys():
                dict_inter[graphe[i][0]] = [graphe[i][1]]
            elif graphe[i][1] in dict_inter[graphe[i][0]]:
                continue
            else:
                dict_inter[graphe[i][0]].append(graphe[i][1])
                
               
        # élimine les arêtes multiples
        
        for key0 in dict_inter.keys():
            dict_inter[key0] = list(set(dict_inter[key0]))
            set_elt = set(dict_inter[key0]).intersection(set(dict_inter.keys()))
            for elt in set_elt:
                if key0 in dict_inter[elt]:
                    dict_inter[elt].remove(key0)
        
        self.writeCleanFile(arg0, dict_inter)

        

        
if __name__ == "__main__":
    
    arg0 = sys.argv[1]
    #arg0 = "/home/benoit/Densite_20/GraphSim1.gr"
    #arg0 = "/home/benoit/Bureau/Stage_M2/Guenocheries/BootClust/Data/HQ08.gr"
     
    CleanGraph(arg0)
            
        

