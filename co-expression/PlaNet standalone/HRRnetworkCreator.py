from Tkinter import *
from tkFileDialog import askopenfilename
from numpy import array, shape, sum, reshape, sqrt, dot, transpose
import time, codecs
class App:
    def __init__(self, master):
        menu = Menu(root)
        root.config(menu=menu)
        self.messenger = Text(root, font=("Courier", 10),wrap=WORD)
        self.printer("Compiler of Pearson Correlation Ranks and Highest Reciprocal Rank (hrr) networks.\nYou can use this tool to generate hrr based networks from expression data.\nThis process requires two steps:\n1)Calculate ranks from expression data\n2)Calculate hrr network from ranks\nUse buttons on the bottom to commence the analysis, but make sure that the expression data follows the format specified in readme.txt.\n")

        Network = LabelFrame(root, text="2. Calculate HRR based co-expression network.")
        DescNetwork = Label(Network, text="This will generate .hrr file containing edge \nand edge-weight information.")
        DescNetwork.pack()
        buttonNetwork = Button(Network, text="Calculate HRR values", fg="red", command=self.loadRANK)
        buttonNetwork.pack()
        Network.pack(side=BOTTOM,anchor=S)

        HRR = LabelFrame(root, text="1. Generate ranks of correlation from expression matrix.")
        DescHRR = Label(HRR, text="This will generate RANK file which is used\nto generate Highest Reciprocal Rank (HRR) networks.\nThis may take minutes/hours, depending on the size of dataset.")
        DescHRR.pack()
        buttonHRR = Button(HRR, text="Load Expression matrix \nand calculate ranks values", fg="red", command=self.LoadExpMat)
        buttonHRR.pack()
        HRR.pack(side=BOTTOM,anchor=N)

        ybar = Scrollbar(root)
        xbar = Scrollbar(root, orient=HORIZONTAL)
        ybar.pack(side=RIGHT, fill=Y)
        xbar.pack(side=BOTTOM, fill=X)
        self.messenger.pack(side=RIGHT,expand=1,fill=BOTH)      
        self.messenger.config(yscrollcommand=ybar.set)
        self.messenger.config(xscrollcommand=xbar.set)
        ybar.config(command=self.messenger.yview)
        xbar.config(command=self.messenger.xview)
        
    def LoadExpMat(self):
        database = askopenfilename(filetypes=[("Expression matrix",".exp")])
        if database !='':
            self.printer("Checking expression matrix integrity...\n")
            DB = codecs.open(database, mode='r', encoding='ASCII', errors='ignore')
            error,size = 0,0
            nominators,denominators= [],[]
            for line in DB:
                    splitta = line.rstrip().split("\t")
                    if size==0:
                        size = len(splitta)
                    else:
                        if size!=len(splitta):
                            self.printer("Warning! Uneven number of columns found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line)
                    if len(splitta)==size:
                        temp = []
                        for j in range(5, len(splitta)):
                            try:
                                temp.append(float(splitta[j]))
                            except:
                                self.printer("Warning! Non-number character found in line:\n%s.\nExpression matrix corrupt. Aborting!\n" % line)
                                error=1
                                break
                        expRow = array(temp)
                        nomi = expRow-(sum(expRow)/len(expRow))
                        nominators.append(nomi)
                        denominators.append(sqrt(sum(nomi**2)))

            if error==0:
                nominators = array(nominators)
                denominators = array(denominators)
                v = open(database.split(".")[0]+".ranks","w")
                v.close()
                v = open(database.split(".")[0]+".ranks","a")
                self.printer("Database OK.\nCalculating Pearson Correlation Coefficient and ranks.\n") 
                for i in range(len(nominators)):
                    print "Calculated PCC values for probeset:%s out of %s." % (str(i+1), len(nominators))
                    dicto = {}
                    nominator = dot(nominators, nominators[i])
                    denominator = dot(denominators, denominators[i])
                    rValues = nominator/denominator
                    for j in range(len(rValues)):
                        dicto[rValues[j]] = str(j)+"\t"
                    sort = dicto.keys()
                    sort.sort(reverse=True)
                    temp = []
                    numOfProbesets = len(sort)
                    for j in range(1,1000):
                        temp.append(dicto[sort[j]])
                    for j in range(numOfProbesets-1000, numOfProbesets):                    
                        temp.append(dicto[sort[j]])           
                    v.writelines(temp+["\n"])                
                v.close()
                self.printer("Ranks calculated and saved as %s." % database.split(".")[0]+"ranks")
            
    def loadRANK(self):
        rankFileName = askopenfilename(filetypes=[("rank of correlation",".ranks")])
        annotations = codecs.open(rankFileName.split(".")[0]+".EXP", mode='r', encoding='ASCII', errors='ignore')
        famMapDic = {}
        indexa=0
        for line in annotations:
            try:
                split = line.rstrip().split("\t")
                famMapDic[str(indexa)] = [split[0], split[0]+"\t"+split[1]+"\t"+split[2]+"\t"+split[3]+"\t"+split[4]+"\t"]
                indexa += 1
            except:
                pass
        
        rankFile = codecs.open(rankFileName, mode='r', encoding='ASCII', errors='ignore')
        netDic = {}
        indexa = 0
        for line in rankFile:            
            split = line.rstrip().split("\t")
            topX = split[:40]
            netDic[indexa] = [topX,[]]
            indexa+=1
        genes = netDic.keys()
        genes.sort()
        
        indexa = 0
        v = codecs.open("%s.hrr" % rankFileName.split(".")[0] ,mode="w", encoding='ASCII', errors='ignore')
        for j in genes:
            print "Calculated HRR values for probeset:%s out of %s." % (indexa, len(genes))
            currentGeneTop = netDic[int(j)][0]
            temp=famMapDic[str(j)][1]
            indexa+=1
            for k in range(len(currentGeneTop)):
                reverseCurrent = netDic[int(currentGeneTop[k])][0]
                try:
                    test = 0
                    test = reverseCurrent.index(str(j))
                    maxa = max(test,k)
                    temp+=str(currentGeneTop[k])+"+%s\t" % maxa
                except  ValueError:
                    pass
                
            currentGeneBottom = netDic[j][1]
            v.write(temp.rstrip()+"\n")

        v.close()
        self.printer("\nHRR based network calculated and saved as %s" % rankFileName.split(".")[0]+".hrr")

    def printer(self, tekst):   #prints text to message box
        self.messenger.insert(END, tekst)
                             
    def clear(self):    #clears message box
        self.messenger.delete(1.0, END)

root = Tk()
root.title("Compiler of co-expression networks for PlaNet")
app = App(root)
root.mainloop()
