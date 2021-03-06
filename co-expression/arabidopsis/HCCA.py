##This script is written by Marek Mutwil.
##This script contains implementation of HCCA algorithm.
name = "ExpMatAra.HRR"   #name of the network file
step = 3                 #step size parameter for generating NVNs
sizeRange = [40,200]     #desired cluster size


a = open(name,"r").readlines()
scoreDic={} 
curDic={}   
loners=[]   
for i in range(len(a)):   ###this loop processess the network file into 2 dictionaries, curDic and scoreDic. 
    splitted = a[i].split("\t")   
    dicto={}
    connections = []
    for j in range(5,len(splitted)):
        if "+" in splitted[j]:
            splitx = splitted[j].split("+")
            if float(splitx[1])<30:
                dicto[splitx[0]]=1/(float(splitx[1])+1)
                connections.append(splitx[0])
    if len(dicto)!=0:
        scoreDic[str(i)] = dicto
        curDic[str(i)] = connections
    else:
        loners.append(str(i))

clustered=[]

clustets=[]
def clustettes(lista):   ##this recursive function extracts isles of nodes that are smaller than...
    cons=[]
    for j in range(len(lista)):
        cons+=curDic[lista[j]]
    cons = list(set(cons+lista))
    if len(cons)>200:    ##...200 nodes
        return
    elif len(cons)==len(lista): 
        cons.sort()
        if cons not in clustets:
            clustets.append(cons)
        return
    else:
        clustettes(cons)

notLoners=curDic.keys()
for i in range(len(notLoners)): ##this loop calls clustettes function to extract and remove nodes forming small islands
    clustettes([notLoners[i]])



for i in range(len(clustets)): ##this loop removes nodes found by clustettes from the main function
    for j in range(len(clustets[i])):
        del curDic[clustets[i][j]]

class CCA(): ##this is the HCCA class
    def SurroundingStep(self, lista, whole,step):  ##this function generates NVN of n steps
        if step<3:  ##step size is defined to 3
            nvn=lista            
            for j in range(len(lista)):
                nvn+=curDic[lista[j]]
            nvn = list(set(nvn))
            self.SurroundingStep(nvn,whole,step+1)            
        else:
            whole.append(lista)
    
    def Chisel(self, NVN, clusters): ##this function recursively removes nodes from NVN. Only nodes that are connected more to the inside of NVN are retained
        temp=[]
        seta = set(NVN)
        for i in range(len(NVN)):
            connections = curDic[NVN[i]]
            inside = set(NVN)&set(connections)
            outside = (set(connections)-set(inside))
            inScore=0
            outScore=0
            for j in inside:
                inScore+=scoreDic[NVN[i]][j]
            for j in outside:
                outScore+=scoreDic[NVN[i]][j]
            if inScore>outScore:
                temp.append(NVN[i])
        if len(temp)==len(seta):
            clusters.append(temp)
            return
        else:
            self.Chisel(temp,clusters)
        
    def BiggestIsle(self, lista, clusterSet,curSeed): ##sometimes the NVN is split into to islands after chiseling. This function finds the biggest island and keeps it. The smaller island is discarded. 
        temp = []
        for k in range(len(lista)):
            temp += scoreDic[lista[k]].keys()
        nodes = set(temp+lista)&clusterSet
        if len(set(nodes))==len(set(lista)):
            curSeed.append(list(set(nodes)))
            return
        else:
            self.BiggestIsle(list(nodes), clusterSet,curSeed)

    def nonOverlappers(self, clusters): ##This function accepts a list of Stable Putative Clusters and greedily extracts non overlapping clusters with highest modularity.
        rankedClust=[]
        for i in range(len(clusters)):
            inScore=0
            outScore=0
            for j in range(len(clusters[i])):
                connections = set(scoreDic[clusters[i][j]].keys())
                inCons = list(connections&set(clusters[i]))
                outCons = list(connections-set(clusters[i]))
                inScore = 0
                outScore = 0
                for k in range(len(inCons)):
                    inScore+=scoreDic[clusters[i][j]][inCons[k]]
                for k in range(len(outCons)):
                    outScore+=scoreDic[clusters[i][j]][outCons[k]]
            rankedClust.append([outScore/inScore, clusters[i]])
            
        rankedClust.sort()
        BestClust = [rankedClust[0][1]]
        for i in range(len(rankedClust)):
            counter=0
            for j in range(len(BestClust)):
                if len(set(rankedClust[i][1])&set(BestClust[j]))>0:
                    counter+=1
                    break
            if counter==0 and rankedClust[i][0]<1:
                BestClust.append(rankedClust[i][1])
        return BestClust

    def networkEditor(self, clustered): ##This function removes nodes in accepted clusters from the current network. 
        connected=[]
        clusteredNodes=[]
        for i in range(len(clustered)):
            clusteredNodes+=clustered[i]
            for j in range(len(clustered[i])):
                connected+=curDic[clustered[i][j]]
                del curDic[clustered[i][j]]
        connections = list(set(connected)-set(clusteredNodes))
        for i in range(len(connections)):
            curDic[connections[i]] = list(set(curDic[connections[i]])-set(clusteredNodes))

    def __init__(self):  ##This function initiates CCA functions
        save=[]
        notClustered = curDic.keys()
        for i in range(len(notClustered)):
            print "Generating SPC for node " +str(i)+" out of " +str(len(notClustered))
            whole=[]
            clusters=[]
            self.SurroundingStep([notClustered[i]],whole, 0)
            self.Chisel(whole[0],clusters)
            if len(clusters[0])>20:
                checked=[]
                for j in range(len(clusters[0])):
                    if clusters[0][j] not in checked:
                        curSeed=[]
                        self.BiggestIsle([clusters[0][j]], set(clusters[0]),curSeed)
                        checked+=curSeed[0]
                        if sizeRange[1]>=len(curSeed[0])>=sizeRange[0]: ##Here the desired size of clusters is specified
                            save.append(curSeed[0])
                            break
        print "finding non-overlappers"
        newCluster = self.nonOverlappers(save)
        print "Found %s non overlapping SPCs. Making a cluster list" % len(newCluster)
        for i in range(len(newCluster)):
            clustered.append(newCluster[i])
        print "%s clusters are now existing. Started the network edit." % len(clustered)
        self.networkEditor(newCluster)
        print "Finished the edit."



def filler(LeftOvers): ##This function assigns nodes that were not clustered by HCCA to clusters they are having highest connectivity to.
    conScoreMat = [[]]*len(clustered)
    clustera=[]
    print len(LeftOvers)
    if len(LeftOvers)!=0:
        for i in range(len(LeftOvers)):
            for j in range(len(clustered)):
                connections = list(set(scoreDic[LeftOvers[i]].keys())&set(clustered[j]))
                score = 0
                for k in range(len(connections)):
                    score+=scoreDic[LeftOvers[i]][connections[k]]
                conScoreMat[j] = score
            
            topScore = max(conScoreMat)
            if topScore!=0:
                sizeList = []
                for j in range(len(conScoreMat)):
                    if conScoreMat[j] == topScore:
                        sizeList.append([len(clustered[j]), j])
                sizeList.sort()
                clustered[sizeList[0][1]] = clustered[sizeList[0][1]]+[LeftOvers[i]]
                clustera.append(LeftOvers[i])
        LeftOvers = list(set(LeftOvers)-set(clustera))
        return filler(LeftOvers)
    else:
        return

mode="go"
iteration=1
CCA()
while mode=="go":   ##This part initiates HCCA class
    try:
        print "iteration: %s" % iteration
        CCA()
        iteration+=1
    except:
        mode="nogo"
        leftovers = curDic.keys()
        filler(leftovers)
        save=[]
        for i in range(len(clustered)):
            for j in range(len(clustered[i])):
                save.append("%s\t%s\n" % (clustered[i][j], str(i)))
        for i in range(len(clustets)):
            for j in range(len(clustets[i])):
                save.append("%s\ts%s\n" % (clustets[i][j], str(i)))
        for i in range(len(loners)):
            save.append("%s\ts%s\n" % (loners[i], "NA")) 
        v = open(name.split(".")[0]+".HCCA","w")
        v.writelines(save)
        v.close()
