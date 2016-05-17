graphVizPath = "C:/Graphviz/bin/neato"
import os, sys
from pylab import hypergeometric

import subprocess

def stepCollector(network, query, steps,nvnP):
    if steps ==0:
        nvnP.append(query)
    else:
        temp = []
        for i in query:
            temp+=network[i][0]
        query = list(set(temp+query))
        stepCollector(network, query, steps-1,nvnP)

def graphFileCreator(nvnP, GIandAnno, dicto):
    nodes = []
    nodesPajek = "*Vertices %s\n" % len(nvnP)
    header=['graph G {\n\toverlap=scale;\n\tranksep=1.5;\n\tconcentrate=true;\n\tsize="12,10";\n']
    edgesPajek = ""
    edges = []
    for i in range(len(nvnP)):
        nodesPajek+="%s %s\n" % (str(i), GIandAnno[nvnP[i]][0])
        nodes.append('\ta%s [shape=point, URL="nodea%s"];\n' % (nvnP[i], nvnP[i]))        
        inNvn = list(set(dicto[nvnP[i]][0])&set(nvnP))        
        temp = []
        for j in range(len(inNvn)):
            if nvnP[i]>inNvn[j]:
                edge = "\ta%s -- a%s ;\n" % (inNvn[j],nvnP[i])
                if edge not in edges:
                    edges.append(edge)
                    edgesPajek+="%s %s\n" % (GIandAnno[inNvn[j]][0],GIandAnno[nvnP[i]][0])
            else:
                edge = "\ta%s -- a%s ;\n" % (nvnP[i],inNvn[j])
                if edge not in edges:
                    edges.append(edge)
                    edgesPajek+="%s %s\n" % (GIandAnno[nvnP[i]][0],GIandAnno[inNvn[j]][0])
    graph = header+edges+nodes+["}"]

    v = open("network.txt","w")
    v.writelines(graph)
    v.close()
    v = open("network.net","w")
    v.writelines(nodesPajek+edgesPajek)
    v.close()
    v = open("network.sif","w")
    v.writelines(edgesPajek)
    v.close()

def mapmanTermEnrichment(nvn, termList, dicto, mapmanDic):
    nvnMap = []
    for i in nvn:
        nvnMap.append(dicto[i][2])
    mapResults = []
    for i in list(set(nvnMap)):
        nGoodInNvn = nvnMap.count(i)
        nGoodTotal = termList.count(i)
        s = hypergeometric(nGoodTotal, len(termList)-nGoodTotal, len(nvnMap), 1000)
        pval = sum(s>=nGoodInNvn)/1000.
        if pval<=0.05:
            try:
                mapResults.append("<tr><td>"+i+"</td><td>"+mapmanDic[i]+"</td><td>"+str(nGoodInNvn)+"</td><td>"+str(nGoodTotal)+"</td><td>"+str(pval)+"</td></tr>\n")
            except:
                pass
    mapResults.sort()
    return mapResults 

def makeNetwork(query, steps, hrr):
    #open network file
    a = open("dbconf.txt","r").read().rstrip()
    query = query.replace("-","_")
    genomeNetwork = open(a.split(".")[0]+".hrr","r").read().replace("-","_").splitlines()
    GIandAnno = {}
    for i in range(len(genomeNetwork)):
        splitted = genomeNetwork[i].replace(";",".").split("\t")
        #GIandAnno["a"+splitted[0]] = [splitted[1], splitted[2].replace("*"," ")]
        GIandAnno[i] = [splitted[0], splitted[1], splitted[2].replace("*"," ")]

    #read in the network
    dicto = {}
    edgeDict = {}
    termList = []
    for i in range(len(genomeNetwork)):
        splitta = genomeNetwork[i].rstrip().split("\t")
        if query == splitta[0]:
            queryIndex = i
        posiC = []
        posiS = []
        negaC = []
        negaS = []
        for j in splitta[5:]:
            if "+" in j:
                splitP = j.split("+")
                if int(splitP[1])<=hrr:
                    #posiC.append("a"+splitP[0])
                    posiC.append(int(splitP[0]))
                    #temp = ["a"+splitta[0], "a"+splitP[0]]
                    temp = [i, int(splitP[0])]
                    temp.sort()
                    if float(splitP[1]) <= int(hrr)/3:
                        edgeDict[temp[0],temp[1]] = "green"
                    elif float(splitP[1]) <= (int(hrr)/3)*2:
                        edgeDict[temp[0],temp[1]] = "orange"
                    elif float(splitP[1]) <= int(hrr):
                        edgeDict[temp[0],temp[1]] = "red"
            if "!" in j:
                splitN = j.split("!")
                if int(splitN[1])<=hrr:
                    #negaC.append("a"+splitN[0])
                    negaC.append(int(splitN[0]))
                    #temp = ["a"+splitta[0], "a"+splitN[0]]
                    temp = [i, int(splitN[0])]
                    temp.sort()
                    if float(splitN[1]) <= int(hrr)/3:
                        edgeDict[temp[0],temp[1]] = "blue"
                    elif float(splitN[1]) <= (int(hrr)/3)*2:
                        edgeDict[temp[0],temp[1]] = "lightblue"
                    elif float(splitN[1]) <= int(hrr):
                        edgeDict[temp[0],temp[1]] = "gray"
        #dicto["a"+splitta[0]] = [posiC,negaC,splitta[4]]
        dicto[i] = [posiC,negaC,splitta[4]]
        termList.append(splitta[4])
    nvnP = []
    stepCollector(dicto, [queryIndex], steps, nvnP)
    nvnP = nvnP[0]
    graphFileCreator(nvnP, GIandAnno, dicto)

    ##Find nodes correlating with negative hhr rank to nodes found in nvn
    nvnN = []
    for i in nvnP:
        nvnN+=dicto[i][1]

    ##MapMan term enrichment calculation
    mapmanDic= {}
    x = open("MapTerms.txt","r").readlines()
    for i in range(len(x)):
        splitta = x[i].rstrip().split("\t")
        mapmanDic[splitta[0]] =splitta[1]
    mapTermEnrP = mapmanTermEnrichment(nvnP, termList, dicto, mapmanDic)    
    v = open("mapRes.mapman","w")
    v.writelines(['<table border="1">\n<tr><td>MapMan term ID</td><td>Term Description</td><td>Terms in NVN</td><td>Terms total</td><td>p-value</td></tr>\n']+mapTermEnrP+["\n</table>\n"])
    v.close()
    
    ###Call graphviz to calculate layout and draw the network
    value = os.system(graphVizPath+" -Gmodel=subset -Timap network.txt -o layout.map")
    if value==1:
        print "Install Graphviz. If Graphviz is installed, correct the path to Graphviz.\n"
    layout = open("layout.map","r").read().replace("a","").splitlines()
    nodes=[]
    xCor = []
    yCor = []    
    for i in range(1,len(layout)):
        s = layout[i].split()
        xCor.append(float(s[2].split(",")[0]))
        yCor.append(float(s[2].split(",")[1]))
        nodes.append(int(s[1].replace("node","")))                
    xmin,xmax = min(xCor), max(xCor)
    ymin,ymax = min(yCor), max(yCor)
    coordDic = {}
    for i in range(1, len(layout)):
        s = layout[i].split()
        coordinata = s[2].split(",")
        coef = 1000/max(xmax,ymax)
        cx = str(float(coordinata[0])*coef-xmin+20) 
        cy = str(float(coordinata[1])*coef-ymin+20)
        coordDic[int(s[1].replace("node",""))] = [cx,cy, coordinata[2]]

    network = open("network.txt","r").read().replace("a","").splitlines()
    edges = []
    for i in range(len(network)):
        if " -- " in network[i]:
            #try:
            s = network[i].split()
            edges.append([int(s[0]),int(s[2]), edgeDict[int(s[0]),int(s[2])]])
            #except:
               # pass

    nodes = list(set(nodes))
    nodeDic={}
    for i in range(len(nodes)):
        temp=[]
        for j in range(len(edges )):
            if nodes[i] == edges [j][0]:
                temp.append(edges [j][1])
            if nodes[i] == edges [j][1]:               
                temp.append(edges [j][0])
        nodeDic[nodes[i]]= temp

    newHead = open("SVG HEAD NEW.txt","r").read()
    whole= newHead %  (xmin, ymin,xmax-xmin, ymax-ymin)

    ids = '''<g id="node%s" class="node" onmouseover="nodes('%s','%s');" onmouseout="denodes('%s','%s');"><title>%s</title>\n'''
    node = '<circle cx="%s" cy="%s" r="6" style="fill:%s" id="%s"/>\n'
    text = '<text text-anchor="middle" x="%s" y="%s" style="font-family:Helvetica;font-size:7.00pt;" onmouseover="ShowTooltip(evt)" onmouseout="HideTooltip(evt)">%s\n'
    texto = '<desc><![CDATA[%s]]></desc>\n'
    button = '<text x="320" y="120">click</text>\n<circle cx="310" cy="118" r="6" onclick="finder()"/>'
    line = '<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1.0">\n</line>\n'
    mouseOver = '<set attributeName="fill" from="lightgrey" to="green" begin="mouseover" end="mouseout"/>\n'
    mOverChanger='<set attributeName="fill" from="black" to="red" begin="node%s.mouseover" end="node%s.mouseout"/>\n'    

    for con in edges: 
        whole+=line % (coordDic[con[0]][0], coordDic[con[0]][1], coordDic[con[1]][0], coordDic[con[1]][1], con[2])                       

    wikiPheno = ""
    x = coordDic.keys()

    for i in range(len(x)):       
        conString="node%s|"
        currentNodeCons=""
        for nodeInd in nodeDic[x[i]]:
            currentNodeCons+= conString % nodeInd       
        temp = ids % (x[i], "node"+str(x[i]), currentNodeCons[:len(currentNodeCons)-1], "node"+str(x[i]), currentNodeCons[:len(currentNodeCons)-1], GIandAnno[x[i]][0])
        if x[i] == queryIndex:
            temp += node % (coordDic[x[i]][0],coordDic[x[i]][1], "aqua", GIandAnno[x[i]][0])
        else:            
            temp += node % (coordDic[x[i]][0],coordDic[x[i]][1], "lightgrey", GIandAnno[x[i]][0])
        temp += text % (coordDic[x[i]][0],coordDic[x[i]][1], GIandAnno[x[i]][1])
        temp+= texto % (GIandAnno[x[i]][2])
        whole+=temp+"</text>\n"+"\n"+"</g>"+"\n"

    v = open("network.svg", "w")
    tipboxsvg = open("tooltip svg.txt","r").read()
    v.write(whole+button+"</g>\n%s\n</svg>" % tipboxsvg)
    v.close()
#makeNetwork("Os.7158.1.S1_a_at", 2, 30)
