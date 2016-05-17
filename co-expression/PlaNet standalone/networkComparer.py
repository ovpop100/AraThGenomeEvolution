graphVizPath = "C:/Graphviz/bin/neato"
import os, random, sys

sys.setrecursionlimit(50000)

def stepCollector(network, query, steps,nvnC, genesFamDict, nvnsByFam, famsInNvns,edgeList,tempEdge,tempNvn): ###This function is used to generate vicinity networks.
    if steps ==0:
        inFams = []
        nvnByFam = []
        for i in query:
            try:
                famsInNvns.append(genesFamDict[i])
                nvnByFam.append(genesFamDict[i])
                inFams.append(i)                
            except:
                pass
        nvnC.append(inFams)
        nvnsByFam.append(nvnByFam)
        edgeList.append(tempEdge)
        tempNvn.append(inFams)
    else:
        temp = []
        for i in query:
            temp+=network[i]
            for j in network[i]:
                try:
                    tempEdgeForStep = [genesFamDict[i],genesFamDict[j]]
                    tempEdgeForStep.sort()
                    if tempEdgeForStep not in tempEdge:
                        tempEdge.append(tempEdgeForStep)
                except:
                    pass
        query = list(set(temp+query))
        stepCollector(network, query, steps-1,nvnC,genesFamDict, nvnsByFam, famsInNvns,edgeList, tempEdge,tempNvn)

def coexpressionGrouper(nvnC,coexpG): ### This function is used to detect and group overlapping vicinity networks.
    print len(coexpG)
    if len(nvnC) == 1:
        score = 0
        for i in range(len(coexpG)):
            if len(set(nvnC[0])&set(coexpG[i]))>1:
                hit = list(set(nvnC[0]+coexpG[i]))
                score=1
                break
        if score==0:
            coexpG.append(nvnC[0])
        else:
            coexpG.remove(coexpG[i])
            coexpG.append(nvnC[0])
    else:
        inGroup = nvnC[0]
        score = 0
        for i in range(1, len(nvnC)):
            if len(set(nvnC[0])&set(nvnC[i]))>0:
                inGroup+=nvnC[i]
                score+=1
                break
        if score==0:
            coexpG.append(list(set(inGroup)))
            coexpressionGrouper(nvnC[1:],coexpG)
        if score>0:
            nvnC.remove(nvnC[i])
            nvnC = nvnC[1:]+[list(set(inGroup))]
            coexpressionGrouper(nvnC ,coexpG)
        pass

def start(step, hrr, networks, query, queryFam):
    ####Read in the networks
    genesFamDict,inFamDict,agiDict = {},{},{}
    famContent,famsInNvns, nvnsByFam, nvnC,edgeList,  = [],[],[],[],[] 
    tableSave = "START OF QUERIES:\n"
    currentSize = 0
    for i in range(len(networks)):
        network,tempFamDict = {},{}
        tempContent = []
        a = open(networks[i],"r").readlines()
        for j in range(len(a)):
            splitta = a[j].rstrip().split("\t")
            temp = []
            agiDict[j+currentSize]= [splitta[1],i, splitta[0]]
            if splitta[3]!="No hit":
                if query == splitta[0]:
                    queryIndex = j+currentSize
                tempContent.append(splitta[3])
                genesFamDict[j+currentSize] = splitta[3]
                try:
                    inFamDict[splitta[3]]+=[j+currentSize]
                except:
                    inFamDict[splitta[3]]=[j+currentSize]
                try:
                    tempFamDict[splitta[3]]+=[j+currentSize]
                except:
                    tempFamDict[splitta[3]]=[j+currentSize]
            for k in range(5,len(splitta)):
                if "+" in splitta[k]:
                    splittb = splitta[k].split("+")
                    if int(splittb[1])<=hrr:
                        temp.append(int(splittb[0])+currentSize)

            network[j+currentSize]= temp
        famContent.append(tempContent)
        famOfQuery = tempFamDict[queryFam]
        for j in famOfQuery:
            tableSave+="%s\t" % agiDict[j][2]
            ####Generate vicinity networks for each member of query's family, in specified networks/organisms
            tempNvn = []
            stepCollector(network, [j], step,nvnC, genesFamDict, nvnsByFam, famsInNvns, edgeList, [],tempNvn)
        currentSize += len(network)

    ####Prepare family/gene ID table for the second step of network comparer (ancestral network)
    tableSave = [tableSave.rstrip()+"\n"]
    famsInNvns = list(set(famsInNvns))
    table = []
    famEnrTable = []
    for i in range(len(famsInNvns)):
        temp = []
        count = 0
        for j in range(len(nvnC)):
            hits = list(set(inFamDict[famsInNvns[i]])&set(nvnC[j]))
            temp.append(hits)
            if len(hits)>0:
                count+=1
        if count>1:
            table.append([count, famsInNvns[i], temp])
    table.sort(reverse=True)

    famEnrTable = '<table border="1">\n<tr><td>Family:</td><td>Description</td><td>Times found in analyzed vicinity networks</td></tr>\n'
    pfamDescDic = {}
    pfamDesc = open("pfamDesc.txt","r").readlines()
    for i in pfamDesc:
        splitta = i.rstrip().split("\t")
        pfamDescDic[splitta[0].replace("-","_")] = splitta[1]
    tableSave += ["START OF TABLE:\n"]
    for i in table:
        try:
            famEnrTable+="<tr><td>"+i[1]+"</td><td>"+pfamDescDic[i[1].replace("-","_")]+"</td><td>"+str(i[0])+"</td></tr>\n"
        except:
            famEnrTable+="<tr><td>"+i[1]+"</td><td> - </td><td>"+str(i[0])+"</td></tr>\n"
        name = i[1]+"\t"
        for j in i[2]:
            column = ""
            if j==[]:
                name+="!"
            else:
                for k in j:
                    column+=agiDict[k][0]+"*"
            name+=column+"\t"
        tableSave.append(name.rstrip()+"\n")
    famEnrTable+="</table>\n"
    ####Prepare edges for drawing ancestral network for the second step of network comparer (ancestral network)
    tableSave.append("START OF EDGES:\n")
    for i in range(len(edgeList)):
        temp = ""
        for j in edgeList[i]:
            temp+=j[0]+"*"+j[1]+"!"
        temp+="\t"
        for j in list(set(nvnsByFam[i])):
            #print nvnsByFam[i]
            temp+=j+"!"
        tableSave.append(temp.rstrip()+"\n")

    ####Detect co-expression groups
    famOfQuery = inFamDict[queryFam]
    coexpG, groupsByProbeset, groupsByFam, groupsForTable = [],[],[],[]
    coexpressionGrouper(nvnC, coexpG)
    for i in range(len(coexpG)):
        tempa = ""
        if queryIndex in coexpG[i]:
            queryIndexByGroup = i
        hits = list(set(famOfQuery)&set(coexpG[i]))
        groupsByProbeset.append(hits)
        temp = []
        for j in coexpG[i]:
            tempa+="%s " % agiDict[j][2]
            temp.append(genesFamDict[j])
        groupsForTable.append(tempa) 
        groupsByFam.append(temp)

    ####Plot co-expression group similarity network graph
    header='graph G {\n\toverlap=scale;\n\tranksep=1.5;\n\tconcentrate=true;\n'
    nodes = ""
    groupScores = {}    
    for i in range(len(groupsByFam)):
        nodes+='\t%s [shape=point, URL="%s"];\n' % ("Group"+str(i+1),"Group"+str(i+1))
        for j in range(i):
            temp = [i, j]
            temp.sort()
            pval = 0
            obsVal = len(set(groupsByFam[i])&set(groupsByFam[j]))
            if obsVal>1:
                for k in range(1000):
                    sample = random.sample(famContent[agiDict[coexpG[j][0]][1]], len(groupsByFam[j]))
                    if len(set(groupsByFam[i])&set(sample))>= obsVal:
                        pval+=1
                if pval/1000.<=0.05:
                    header += "\tGroup%s -- Group%s [penwidth=%s];\n" % (temp[0]+1,temp[1]+1,obsVal)
                    try:
                        groupScores["Group"+str(temp[0]+1)]+="Group"+str(temp[1]+1)+":("+str(obsVal)+","+str(pval/1000.)[:6]+") "
                    except:
                        groupScores["Group"+str(temp[0]+1)]="Group"+str(temp[1]+1)+":("+str(obsVal)+","+str(pval/1000.)[:6]+") "
                    try:
                        groupScores["Group"+str(temp[1]+1)]+="Group"+str(temp[0]+1)+":("+str(obsVal)+","+str(pval/1000.)[:6]+") "
                    except:
                        groupScores["Group"+str(temp[1]+1)]="Group"+str(temp[0]+1)+":("+str(obsVal)+","+str(pval/1000.)[:6]+") "
                    
    v = open("NetCompNet.txt","w")
    v.write(header+nodes+"\n}")
    v.close()

    value = os.system(graphVizPath+" -Gmodel=subset -Timap NetCompNet.txt -o NetCompNet.map")
    if value==1:
        print "Install Graphviz. If Graphviz is installed, correct the path to Graphviz.\n"        
    layout = open("NetCompNet.map","r").readlines()
    nodes=[]
    xCor = []
    yCor = []    
    for i in range(1,len(layout)):
        s = layout[i].split()
        xCor.append(float(s[2].split(",")[0]))
        yCor.append(float(s[2].split(",")[1]))
        nodes.append(s[1])                
    xmin,xmax = min(xCor), max(xCor)
    ymin,ymax = min(yCor), max(yCor)
    coordDic = {}
    for i in range(1, len(layout)):
        s = layout[i].split()
        coordinata = s[2].split(",")
        coef = 1000/max(xmax,ymax)
        cx = str(float(coordinata[0])*coef-xmin+20)
        cy = str(float(coordinata[1])*coef-ymin+20)
        coordDic[s[1]] = [cx,cy]

    network = open("NetCompNet.txt","r").readlines()

    edgeDict = {}
    nodeDic={}
    for i in range(len(network)):
        if " -- " in network[i]:            
            s = network[i].split()
            edgeDict[s[0],s[2]]= float(s[3].split("=")[1].split("]")[0])
            try:
                nodeDic[s[0]]+=[s[2]]
            except:
                nodeDic[s[0]]=[s[2]]
            try:
                nodeDic[s[2]]+=[s[0]]
            except:
                nodeDic[s[2]]=[s[0]]
    maxEdge = max(map(int, edgeDict.values()))

    newHead = open("SVG HEAD NEW.txt","r").read()
    whole= newHead %  (xmin, ymin,xmax-xmin, ymax-ymin)

    ids = '''<g id="node%s" class="node" onmouseover="nodes('%s','%s');" onmouseout="denodes('%s','%s');"><title>%s</title>\n'''
    node = '<circle cx="%s" cy="%s" r="6" style="fill:%s" id="%s"/>\n'
    text = '<text text-anchor="middle" x="%s" y="%s" style="font-family:Helvetica;font-size:7.00pt;" onmouseover="ShowTooltip(evt)" onmouseout="HideTooltip(evt)">%s\n'
    texto = '<desc><![CDATA[Connected to other groups with (score, p-value): %s]]></desc>\n'
    button = '<text x="320" y="120">click</text>\n<circle cx="310" cy="118" r="6" onclick="finder()"/>'
    line = '<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="lightblue" stroke-width="%s">\n</line>\n'
    mouseOver = '<set attributeName="fill" from="lightgrey" to="green" begin="mouseover" end="mouseout"/>\n'
    mOverChanger='<set attributeName="fill" from="black" to="red" begin="node%s.mouseover" end="node%s.mouseout"/>\n'    

    for con in edgeDict.keys(): 
        whole+=line % (coordDic[con[0]][0], coordDic[con[0]][1], coordDic[con[1]][0], coordDic[con[1]][1], (edgeDict[con]/maxEdge)*10)                      


    x = coordDic.keys()

    for i in range(len(x)):       
        conString="node%s|"
        currentNodeCons=""
        currentNodeVals = ""
        try:
            for nodeInd in nodeDic[x[i]]:
                currentNodeCons+= conString % nodeInd       
        except:
            pass
        temp = ids % (x[i], "node"+x[i], currentNodeCons[:len(currentNodeCons)-1], "node"+x[i], currentNodeCons[:len(currentNodeCons)-1], "")
        if "Group"+str(queryIndexByGroup+1) == x[i]:
            temp += node % (coordDic[x[i]][0],coordDic[x[i]][1], "red", x[i])
        else:
            temp += node % (coordDic[x[i]][0],coordDic[x[i]][1], "lightgreen", x[i])
        temp += text % (coordDic[x[i]][0],coordDic[x[i]][1], x[i])
        try:
            temp+= texto % groupScores[x[i]]
        except:
            pass
        whole+=temp+"</text>\n"+"\n"+"</g>"+"\n"

    v = open("nc.svg", "w")
    tipboxsvg = open("tooltip svg.txt","r").read()
    v.write(whole+button+"</g>\n%s\n</svg>" % tipboxsvg)
    v.close()

    ####Calculate the table displaying similarity scores and p-values of all family members vs. query
    results = {}
    queryNvnbyFam = nvnsByFam[famOfQuery.index(queryIndex)]
    for i in range(len(famOfQuery)):        
        obsVal = len(set(queryNvnbyFam)&set(nvnsByFam[i]))
        pval=0
        for k in range(1000):
            sample = random.sample(famContent[agiDict[famOfQuery[i]][1]] , len(nvnsByFam[i]))
            if len(set(queryNvnbyFam)&set(sample))>= obsVal:
                pval+=1
        results[famOfQuery[i]] = [obsVal, pval/1000.]        

    ####Html formatting of the results
    header = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n<html>\n<head>\n<!-- blarg -->\n</head>\n<body>\n'
    header +='<b>Query: %s belongs to family:%s and could be found %s times in selected networks.</b><br>\n' % (query, queryFam, len(famOfQuery))
    header +='Using step size of: %s and hrr cut-off of: %s, the vicinity networks of %s family was grouped into %s groups.<br><br><br>\n' % (step, hrr, queryFam, len(coexpG))
    header +='<b>Table below shows identity of the probeset IDs constituting each group</b>.<br>\n<table border="1">\n<tr><td>Group:</td><td>Network</td><td>Probesets</td></tr>\n'
    for i in range(len(groupsForTable)):
        header+="<tr><td>Group%s</td><td>%s</td><td>%s</td></tr>\n" % (str(i+1), networks[agiDict[groupsByProbeset[i][0]][1]], groupsForTable[i])    
    header += '</table>\n<br><br><b>Comparison of the groups in respect to the gene family content is shown in the network below</b>.<br>Any two groups that have significantly (p<0.05) similar gene family content are connected. The thickness of connecting edge is proportional to the score (number of gene families in common).<br>\nGroup containing the query gene is colored red.<br>\n'
    header += '<embed src="nc.svg" width="1200" height="1200" type="image/svg+xml" pluginspage="http://www.adobe.com/svg/viewer/install/" />\n<br>'

    #header += '<iframe src="nc.svg" width="1500" height="1500">\n</iframe><br>\n'
    header += '<b>The table below shows score of your query gene to all genes found in the family: %s</b>' % genesFamDict[queryIndex]
    header += '<table border="1">\n<tr><td>Group:</td><td>Network:</td><td>Probeset ID</td><td>Transcript ID</td><td>Score</td><td>p-value</td></tr>\n'
    tableSave.append("START OF VS. QUERY:\n")

    for i in range(len(groupsByProbeset)):
        hiddenTable = ""
        table="<tr><td>Group: %s</td><td>%s</td><td>" % (str(i+1), networks[agiDict[groupsByProbeset[i][0]][1]])
        agi = "<td>"
        score="<td>"
        pval="<td>"
        for j in groupsByProbeset[i]:
            table+=agiDict[j][2]+"<br>"
            agi+=agiDict[j][0]+"<br>"
            score+=str(results[j][0])+"<br>"
            pval+=str(results[j][1])+"<br>"
            hiddenTable +=("Group: %s\t" % str(i+1))+agiDict[j][2]+"\t"+agiDict[j][0]+"\t"+str(results[j][0])+"\t"+str(results[j][1])+"\n"
        table+="</td>"+agi+"</td>"+score+"</td>"+pval+"</td></tr>\n"
        header+=table
        tableSave.append(hiddenTable)
    header+= "</table>\n"
    header+= "<br><b>The table below shows which pfam families are enriched in vicinity networks of %s family in analyzed species.<br>" % queryFam
    header+=famEnrTable
    header+="<b>The analysis is done. Now proceed to the second step in NetworkComparer pipeline by navigating to Analyze/NetworkComparer step II.<br>\n"
    v = open("NetComp.html","w")
    v.writelines([header]+["<!--\n"]+tableSave+["--!>\n"])
    v.close()
    
#start(2, 30, ["ExpMatAra.hrr","ExpMatBar.hrr"], "253428_at", "Cellulose_synt")


