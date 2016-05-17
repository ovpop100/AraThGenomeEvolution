graphVizPath = "C:/Graphviz/bin/neato"
import os

def start(selected):
    pfamDescDic = {}
    pfamDesc = open("pfamDesc.txt","r").readlines()
    for i in pfamDesc:
        splitta = i.rstrip().split("\t")
        pfamDescDic[splitta[0]] = splitta[1]
    
    a = open("NetComp.html","r").read().replace("-","_").splitlines()   
    table,edges = [],[]
    tabla,edga = 0,0
    for i in range(len(a)):
        if "START OF TABLE" in a[i]:
            tabla = 1
        if "START OF EDGES" in a[i]:
            tabla,edga = 0,1
        if "START OF VS. QUERY" in a[i]:
            break
        if "START OF QUERIES:" in a[i]:
            queries = a[i+1].rstrip().split("\t")
        if tabla ==1:
            table.append(a[i])
        if edga ==1:
            edges.append(a[i])

    tableSave, selectedEdges,selectedNodes,allEdges,allNodes  = [],[],[],[],[]
    #tableHeader="<tr><td>pfam family</td><td>Description</td></tr>\n<tr><td></td>"
    tableHeader = "<tr><td>pfam family</td><td>Description</td>"
    for i in selected:
        splitta = edges[i+1].split()
        tableHeader += "<td>%s</td>" % queries[i]
        if len(splitta)==2:
            selectedEdges.append(edges[i+1].split()[0].split("!"))
            selectedNodes.append(edges[i+1].split()[1].split("!"))        
            allEdges+= edges[i+1].split()[0].split("!")
            allNodes+= edges[i+1].split()[1].split("!")
        if len(splitta)==1:
            selectedEdges.append([])
            selectedNodes.append(edges[i+1].split()[0].split("!"))        
            allEdges+= []
            allNodes+= edges[i+1].split()[0].split("!")            

    tableSave.append(tableHeader+"</tr>")
    for i in range(1, len(table)):
        splitta = table[i].rstrip().split()
        try:
            name = "<tr><td>"+splitta[0]+"</td><td>"+pfamDescDic[splitta[0]]+"</td>"
        except:
            name = "<tr><td>"+splitta[0]+"</td><td>-</td>"
        count =0
        for j in selected:
            if splitta[j+1]!="!":
                count+=1
            name+="<td>"+splitta[j+1].replace("*","<br>").replace("!","-")+"</td>"
        if count>1:
            tableSave.append(name.rstrip()+"</tr>")

    allEdges=list(set(allEdges))
    allEdges.remove("")
    allNodes=list(set(allNodes))
    allNodes.remove("")
    edgeCount = {}
    for i in allEdges:
        count = 0
        for j in selectedEdges:
            if i in j:
                count+=1
        if count>1:
            edgeCount[i] =count

    nodeCount = {}           
    for i in allNodes:
        count = 0
        for j in selectedNodes:
            if i in j:
                count+=1
        if count>1:
            nodeCount[i] = count   

    divisions = list(set([round(len(selected),0), round(len(selected)*(3/4.),0), round(len(selected)/2.,0)]))
    divisions.sort(reverse=True)
    ranges = []
    for i in divisions:
        if i>=2:
            ranges.append(i)
    ranges = [ranges, ["green","orange","red"]]
    rangeTextColor, rangeTextNumber= "",""
    for i in range(len(ranges[0])):
        rangeTextColor += str(ranges[1][i])+", "
        rangeTextNumber += str(ranges[0][i])+", "
    print rangeTextColor, rangeTextNumber
    header=['graph G {\n\toverlap=scale;\n\tranksep=1.5;\n\tconcentrate=true;\n']
    graphNodes, graphEdges = [],[]
    for i in nodeCount.keys():
        if nodeCount[i]>=min(divisions):
            graphNodes.append('\ta%s [shape=point, URL="a%s"];\n' % (i, i))
    for i in edgeCount.keys():
        if edgeCount[i]>=min(divisions):
            edgeComponent = i.split("*")
            graphEdges.append('\ta%s -- a%s ;\n' % (edgeComponent[0], edgeComponent[1]))  

    v = open("ancNet.txt","w")
    v.writelines(header+graphEdges+graphNodes+["}"])
    v.close()
    
    value = os.system(graphVizPath+" -Gmodel=subset -Timap ancNet.txt -o ancNet.map")
    if value==1:
        print "Install Graphviz. If Graphviz is installed, correct the path to Graphviz.\n"
    layout = open("ancNet.map","r").readlines()
    nodes=[]
    xCor = []
    yCor = []    
    for i in range(1,len(layout)):
        s = layout[i].split()
        xCor.append(float(s[2].split(",")[0]))
        yCor.append(float(s[2].split(",")[1]))
        nodes.append(s[1][1:])                
    xmin,xmax = min(xCor), max(xCor)
    ymin,ymax = min(yCor), max(yCor)
    coordDic = {}
    for i in range(1, len(layout)):
        s = layout[i].split()
        coordinata = s[2].split(",")
        coef = 1000/max(xmax,ymax)
        cx = str(float(coordinata[0])*coef-xmin+20)
        cy = str(float(coordinata[1])*coef-ymin+20)
        coordDic[s[1][1:]] = [cx,cy, coordinata[2]]


    network = open("ancNet.txt","r").readlines()
    edges = []
    edgeDict = {}
    for i in range(len(network)):
        if " -- " in network[i]:            
            s = network[i].split()
            for j in ranges[0]:
                if edgeCount[s[0][1:]+"*"+s[2][1:]]>=j:
                    indexa = ranges[0].index(j)
                    edges.append([s[0][1:],s[2][1:], ranges[1][indexa]])
                    break


    nodes = list(set(nodes))
    nodeDic={}
    nodeColor = {}
    for i in range(len(nodes)):
        temp=[]
        for j in range(len(edges )):
            if nodes[i] == edges [j][0]:
                temp.append(edges [j][1])
            if nodes[i] == edges [j][1]:               
                temp.append(edges [j][0])
        nodeDic[nodes[i]]= temp       
        for j in ranges[0]:
            if nodeCount[nodes[i]]>=j:
                indexa = ranges[0].index(j)
                nodeColor[nodes[i]] = ranges[1][indexa]
                break

    newHead = open("SVG HEAD NEW.txt","r").read()
    whole= newHead %  (xmin, ymin,xmax-xmin, ymax-ymin)

    ids = '''<g id="node%s" class="node" onmouseover="nodes('%s','%s');" onmouseout="denodes('%s','%s');"><title>%s</title>\n'''
    node = '<circle cx="%s" cy="%s" r="6" style="fill:%s" id="%s"/>\n'
    text = '<text text-anchor="middle" x="%s" y="%s" style="font-family:Helvetica;font-size:7.00pt;" onmouseover="ShowTooltip(evt)" onmouseout="HideTooltip(evt)">%s\n'
##    texto = '<desc><![CDATA[%s]]></desc>\n'
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
        temp = ids % (x[i], "node"+x[i], currentNodeCons[:len(currentNodeCons)-1], "node"+x[i], currentNodeCons[:len(currentNodeCons)-1], "")
        temp += node % (coordDic[x[i]][0],coordDic[x[i]][1], nodeColor[x[i]], x[i])
        temp += text % (coordDic[x[i]][0],coordDic[x[i]][1], x[i])
##        temp+= texto % ("")
        whole+=temp+"</text>\n"+"\n"+"</g>"+"\n"

    v = open("ancNet.svg", "w")
    tipboxsvg = open("tooltip svg.txt","r").read()
    v.write(whole+button+"</g>\n%s\n</svg>" % tipboxsvg)
    v.close()

    header = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n<html>\n<head>\n<!-- NetworkComparer -->\n</head>\n<body>\n'
    header += '<b>The network below shows which gene families are over-represented in vicinity network of selected probesets.</b>\n<br>Nodes and edges colored '+rangeTextColor[:len(rangeTextColor)-2]+' represent gene families and co-expression relationships between them, that are present in '+rangeTextNumber[:len(rangeTextNumber)-2]+' of the selected networks.<br>\n'
    header += '<embed src="ancNet.svg" width="1200" height="1200" type="image/svg+xml" pluginspage="http://www.adobe.com/svg/viewer/install/" />\n<br>'
    #header += '<iframe src="ancNet.svg" width="1500" height="1500">\n</iframe>\n'
    header += '<b>The table below displays the probeset identity of enriched families in selected queries.</b><br>\n'
    header += '<table border="1">\n'
    v = open("ancNet.html","w")
    v.writelines([header]+tableSave)
    v.close()

#start([16, 11, 30])
