from Tkinter import Tk, Menu, Text, Button, Radiobutton, Checkbutton, Listbox, StringVar, IntVar, Label, Scrollbar, Toplevel, LabelFrame, Entry, NONE, END, SUNKEN, HORIZONTAL, WORD, BOTTOM, TOP, LEFT, RIGHT, NW, E, N, W, S, Y, X, BOTH, MULTIPLE
from time import time
from tkFileDialog import askopenfilename
from pylab import average, show
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
import shutil,os ,sys, time, codecs, matplotlib, numpy, networkViewer, networkComparer, ancestralNetwork
os.environ['HOME'] = '/expProfile'


class App:
######################################################################## 
####################################Interface for PlaNet: 
    def __init__(self, master):
        menu = Menu(root)
        root.config(menu=menu)
        filemenu = Menu(menu)
        analyzemenu = Menu(menu)
        PostAnalysis = Menu(menu)
        Settings = Menu(menu)
        About = Menu(menu)
        
        menu.add_cascade(label="Start", menu=filemenu)
        filemenu.add_command(label="Select genes for analysis", command=self.SearchGene)       
        menu.add_cascade(label="Analyze", menu=analyzemenu)
        analyzemenu.add_command(label="Display probeset specific report", command=self.nvnReport)        
        analyzemenu.add_command(label="NetworkComparer Step I: Find similar co-expression networks", command=self.netComparer)
        analyzemenu.add_command(label="NetworkComparer Step II: Extract ancestral network", command=self.ancestranNetPipe)
        menu.add_cascade(label="Change organism/database", command=self.UseDatabaseFile)
        menu.add_cascade(label="Help", menu=About)
        About.add_command(label="About GeneCAT", command=self.About)
        
        self.messenger = Text(root, font=("Courier", 10), wrap=NONE)
        self.printer("Stand alone version of PlaNet\n")
        button = Button(root, text="Clear", fg="red", command=self.clear)
        self.LabVar = StringVar()
        DBLabel = Label(root, relief=SUNKEN, textvariable=self.LabVar)
        
        ybar = Scrollbar(root)
        xbar = Scrollbar(root, orient=HORIZONTAL)
        DBLabel.pack(side=BOTTOM, anchor=NW)
        button.pack(side=BOTTOM, anchor=E)
        ybar.pack(side=RIGHT, fill=Y)
        xbar.pack(side=BOTTOM, fill=X)
        self.messenger.pack(expand=1,fill=BOTH)        
        self.messenger.config(yscrollcommand=ybar.set)
        self.messenger.config(xscrollcommand=xbar.set)
        ybar.config(command=self.messenger.yview)
        xbar.config(command=self.messenger.xview)
        
        config = open('dbconf.txt', 'r').readlines()
        self.queries = []
        try:
            self.currentDB = codecs.open('%s' % (config[0].rstrip()), mode='r', encoding='ASCII', errors='ignore').readlines()
            self.printer("Currently using database: %s consisting of %d probesets and %d chips.\n" % (config[0].rstrip(), len(self.currentDB), (len(self.currentDB[0].split())-3)))
            self.LabVar.set("Using database: %s, Columns=%d, rows=%d." %(config[0].rstrip(), len(self.currentDB[0].split()), len(self.currentDB)))
            nominator = []
            denominator = []
            self.annoDict = {}
            for i in range(len(self.currentDB)):
                temp = []
                splitta = self.currentDB[i].rstrip().split("\t")
                self.annoDict[i] = splitta[0]+"\t"+splitta[1]+"\t"+splitta[2].replace("*"," ")+"\t"+splitta[3]+"\t"
                for j in range(5,len(splitta)):
                    temp.append(float(splitta[j]))                    
                expVector = numpy.array(temp)
                nomi = expVector-(numpy.sum(expVector)/len(expVector))
                denomi = numpy.sqrt(numpy.sum(nomi**2))
                nominator.append(nomi)
                denominator.append(denomi)
            self.nominator = numpy.array(nominator)
            self.denominator = numpy.array(denominator)
        except IOError:
            self.printer("Could not open current database.\n")
            
######################################################################## 
####################################Start menu items: 
    def SearchGene(self): #interface for "Select genes for analysis" 
        SearchFrame = Toplevel()
        SearchFrame.title("Data entry. Enter gene identifier, probesets or keywords.\n")
        scrollBar = Scrollbar(SearchFrame)        
        self.inputFrame = Text(SearchFrame, font="courier", width=30, wrap=WORD)
        self.inputFrame.pack(side=LEFT, fill=BOTH)
        scrollBar.pack(side=LEFT, fill=Y)
        self.inputFrame.config(yscrollcommand=scrollBar.set)
        scrollBar.config(command=self.inputFrame.yview)
        goButton = Button(SearchFrame, text="Find", fg="red", command=self.SearchGenePipe)
        goButton.pack(side=BOTTOM)
####################################Searches the current database with user specified keywords
    def SearchGenePipe(self): 
        if self.inputFrame.get(1.0, END)!="\n":
            queries = self.inputFrame.get(1.0, END).lower().split()
            foundGenes = []
            for i in queries:
                infos = ""
                for j in self.currentDB:
                    if i in j.lower():
                        splitta = j.split("\t")
                        foundGenes.append(j)
                        infos+=splitta[0]+"\t"+splitta[1]+"\t"+splitta[2].replace("*"," ")+"\n"
                self.printer("Probesets found for term: "+i+"\n"+infos+"\n")
            if len(infos)=="":
                self.printer("No hits found.\n")
            else:
                self.queries = foundGenes

######################################################################## 
####################################Data mining menu items
                
####################"Display probeset specific report"
    def nvnReport(self): #interface for "Display probeset specific report" 
        top = Toplevel()
        top.title("Co-expression analysis")
        self.genelist=[]
        for string in self.queries:
            self.genelist.append(string.split())            
        self.listbox = Listbox(top, width=40, height=30,exportselection=0)
        for gene in self.genelist:
            self.listbox.insert(END, gene[0]+'   '+gene[1])
            
        DescriptionLabel = LabelFrame(top, text="Info")
        Description = Label(DescriptionLabel, text="Select gene of interest from the list to the left.\nThis tool will generate result.html file containing a page similar to\ngene specific pages in PlaNet.")
        DescriptionLabel.grid(row=0, column=2)
        
        ParametersLabel = LabelFrame(top, text="Parameters")
        Steps = Label(ParametersLabel, text="Number of steps")
        ParametersLabel.grid(row=1, column=2)
        Hrr = Label(ParametersLabel, text="HRR cut-off.")
        
        self.paraSteps = Entry(ParametersLabel)
        self.paraSteps.grid(row=1)
        self.paraSteps.insert(END, 2)
        self.paraHrr = Entry(ParametersLabel)
        self.paraHrr.grid(row=3)
        self.paraHrr.insert(END, 30)      
        Description.grid()
        Steps.grid(row=0)
        Hrr.grid(row=2)
        scrollbar = Scrollbar(top)
        scrollbar.grid(row=0, column=1, rowspan=5, sticky=S+N)   
        self.listbox.grid(row=0, column=0, rowspan=5)
        scrollbar.config(command=self.listbox.yview)
        button = Button(top, text="Calculate!", fg="red", font=("Courier", 22), command=self.nvnReportPipe)
        button.grid(row=6, column=0)

    def nvnReportPipe(self): #Executes functions that generate result.html file
        try:
            steps = int(self.paraSteps.get())
            hrr = int(self.paraHrr.get())
        except:
            self.printer("You need to enter integer values for HRR and step parameters.\n")
        if self.listbox.curselection()!=():
            query = self.genelist[int(self.listbox.curselection()[0])][0]

            #### Plot expression profile function
            try:
                plotDatabase = codecs.open(open('dbconf.txt', 'r').read().rstrip().split(".")[0]+".plt",mode='r', encoding='ASCII', errors='ignore').readlines()
                ticka = [""]+plotDatabase[0].replace(",","\n").lstrip().rstrip().split("\t")
                for i in plotDatabase:
                    if query in i:
                        query = i                
                data = query.split()[1:]
                temp = []
                for i in range(len(data)):
                    temp.append([ map(float,data[i].replace("-","\t").rstrip().split()), average(map(float,data[i].replace("-","\t").rstrip().split()))] )            
                fig = plt.figure(figsize=(12,7))
                ax = fig.add_subplot(111)
                plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.3)
                ax.set_ylabel("Signal value")
                ax.set_title(query.split()[0])
                ax.grid(True)
                plt.xticks(range(len(ticka)+1), ticka, rotation=90, fontsize="small", horizontalalignment="center")
                ax.plot([0],[0])
                crossX = []
                crossY = []
                for i in range(len(temp)):
                    ax.plot([i+1]*len(temp[i][0]),temp[i][0],"g.")
                    crossX.append([i+1])
                    crossY.append(temp[i][1])
                ax.plot(crossX, crossY, "-ro")
                ax.plot([i+2],[0])
                canvas = FigureCanvasAgg(fig)
                canvas.print_figure("profile.png")
                plt.clf()
            except:
                self.printer("Failed to generate an expression profile of your gene of interes.\nThe expression matrix used for plotting of expression profiles must be present and named "+open('dbconf.txt', 'r').read().rstrip().split(".")[0]+".plt!")  
            ###Call network creator
            try:
                networkViewer.makeNetwork(query.split()[0], steps, hrr)
            except:
                self.printer("Failed to generate an co-expression network of your gene of interes.\nThe HRR network file used must be present named "+open('dbconf.txt', 'r').read().rstrip().split(".")[0]+".hrr!")  
            ### Calculate PCC of a gene to all genes in database
            try:
                query = self.queries[int(self.listbox.curselection()[0])].split("\t")
                expVector = map(float, query[5:])
                expVector= numpy.array(expVector)
                nomi = expVector-(numpy.sum(expVector)/len(expVector))
                denomi = numpy.sqrt(numpy.sum(nomi**2))
                rValues = numpy.dot(self.nominator, nomi)/numpy.dot(self.denominator,denomi)
                displayList = []
                for i in range(len(rValues)):
                    displayList.append([rValues[i], self.annoDict[i]])
                displayList.sort(reverse=True)
            except:
                displayList = []
                self.printer("Failed to calculate Pearson correlation co-efficient list.\n")                    

            ###Create html document with results
            header = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n<html>\n<head>\n<!-- blarg -->\n</head>\n<body>\n'
            try:
                header +='<big><b>Summary page for: %s</b></big>\n<br><b>Expression profile:</b>\n<IMG SRC="profile.png"><br>\n' % (displayList[0][1].split()[0]+"\t"+displayList[0][1].split()[1])
            except:
                pass
            header +='<b>HRR based co-expression network:</b>\nGreen, organge and red edges represent HRR values of %s, %s and %s, respectively.</b><br>\n' % (int(hrr)/3, (int(hrr)/3)*2, hrr) 
            header += '<embed src="network.svg" width="1200" height="1200" type="image/svg+xml" pluginspage="http://www.adobe.com/svg/viewer/install/" />\n<br>'
            #header += '<iframe src="network.svg" width="1500" height="1500">\n</iframe>\n'
            header += "<br><b>MapMan ontology analysis of the above network:</b>\n"
            try:
                header += open("mapRes.mapman","r").read()
            except:
                pass
            header+="\n<br><b>Pearson correlation based co-expression analysis:</b>\n<pre>"
            
            v = open("result.html","w")
            v.close()
            v = open("result.html","a")
            v.write(header)
            for i in range(len(self.annoDict)):
                v.write(str(displayList[i][0])[:6]+"\t"+displayList[i][1]+"\n")
            v.write("</pre>")
            v.close()
            self.printer("Probeset specific result calculated and available in result.html file.\n")
            
####################"NetworkComparer step I:"
    def netComparer(self): #Interface
        top = Toplevel()
        top.title("NetworkComparer")
        
        self.genelistNC=[]
        for string in self.queries:
            if "No hit" not in string:
                self.genelistNC.append(string.split("\t"))
        DescriptionLabel = LabelFrame(top, text="Select the gene of interest and specify parameters.")
        Description = Label(DescriptionLabel, text="Select gene of interest from the list below.\nThis tool will generate netComp.html file containing a page similar to\nfirst step of NetworkComparer.")
        DescriptionLabel.grid(row=0, column=0)
        Description.grid(row=0)        
        ParametersLabel = LabelFrame(DescriptionLabel, text="Parameters")
        Steps = Label(ParametersLabel, text="Number of steps")
        ParametersLabel.grid(row=1, column=2)
        Hrr = Label(ParametersLabel, text="HRR cut-off.")
        Steps.grid(row=0)
        Hrr.grid(row=2)        
        self.paraSteps = Entry(ParametersLabel)
        self.paraSteps.grid(row=1)
        self.paraSteps.insert(END, 2)
        self.paraHrr = Entry(ParametersLabel)
        self.paraHrr.grid(row=3)
        self.paraHrr.insert(END, 30)        
        self.listbox = Listbox(DescriptionLabel, width=40, height=30,exportselection=0)
        probeName = []
        for gene in self.genelistNC:
            self.listbox.insert(END, gene[0]+'   '+gene[1]+'   '+gene[3])
        self.listbox.grid(row=1, column=0, rowspan=5, sticky=N+S+E+W)        
        scrollbarG = Scrollbar(DescriptionLabel)
        scrollbarG.grid(row=1, column=1, rowspan=5, sticky=S+N)
        scrollbarG.config(command=self.listbox.yview)
        
        NetworkLabel = LabelFrame(top, text="Select the networks/species you want to compare.")
        Network = Label(NetworkLabel, text="Available networks are displayed in the list below. Check/uncheck networks of interest.")
        NetworkLabel.grid(row=0, column=1)
        Network.grid(row=0)
        self.networkBox = Listbox(NetworkLabel, width=40, height=30, selectmode=MULTIPLE, exportselection=0)
        self.networkList = []
        for i in os.listdir("."):
            if ".hrr" in i:
                self.networkBox.insert(END, i)
                self.networkList.append(i)
        self.networkBox.grid(row=1, column=0, rowspan=5, sticky=N+S+E+W)
        scrollbarN = Scrollbar(NetworkLabel)
        scrollbarN.grid(row=1, column=1, rowspan=5, sticky=S+N)
        scrollbarN.config(command=self.networkBox.yview)
        button = Button(top, text="Calculate!", fg="red", font=("Courier", 22), command=self.netComparerPipe)
        button.grid(row=1, column=0, columnspan=5, sticky=E+W)      

    def netComparerPipe(self): #Passess selected query and networks to networkComparer.start function
        error = 0
        try:
            steps = int(self.paraSteps.get())
        except:
            error+=1
            self.printer("You must enter an integer value bigger than 1 for this parameter.\n")
        try:
            hrr = int(self.paraHrr.get())
        except:
            error+=1
            self.printer("You must enter an integer value bigger than 1 for this parameter.\n")
        netList = []
        for i in self.networkBox.curselection():
            netList.append(self.networkList[int(i)])
        if netList==[]:
            error+=1
            self.printer("You must at least select one network for this analysis.\n")
        if error==0:
            #try:
                networkComparer.start(steps, hrr, netList, self.genelistNC[int(self.listbox.curselection()[0])][0],self.genelistNC[int(self.listbox.curselection()[0])][3])
                self.printer("Calculation done and results are available in NetComp.html file.\n")
            #except:
                #self.printer("Something went wrong. Make sure Graphviz is installed.\n")
                
####################"NetworkComparer step II:"      
    def ancestranNetPipe(self): #Interface
        top = Toplevel()
        top.title("NetworkComparer")
        a = open("NetComp.html","r").readlines()
        self.queriesAC = []
        quera = 0
        for i in range(len(a)):
            if "START OF QUERIES" in a[i]:
                self.queryOrder = a[i+1].rstrip().split()
            if "START OF VS. QUERY" in a[i]:
                quera=1                
            if quera==1:
                self.queriesAC.append(a[i].rstrip().split("\t"))
        self.queriesAC = self.queriesAC[1:len(self.queriesAC)-1]
        DescriptionLabel = LabelFrame(top, text="Select the gene of interest and specify parameters.")
        Description = Label(DescriptionLabel, text="Select gene of interest from the list below.\nThis tool will generate netComp.html file containing a page similar to\nfirst step of NetworkComparer.")
        DescriptionLabel.grid(row=0, column=0)
        Description.grid(row=0)        
        self.listbox = Listbox(DescriptionLabel, width=40, height=30,exportselection=0,selectmode=MULTIPLE)
        for gene in self.queriesAC:
            self.listbox.insert(END, gene[0]+'   '+gene[1]+'   '+gene[2]+'   '+gene[3]+'   '+gene[4])
        self.listbox.grid(row=1, column=0, rowspan=5, sticky=N+S+E+W)        
        scrollbarG = Scrollbar(DescriptionLabel)
        scrollbarG.grid(row=1, column=1, rowspan=5, sticky=S+N)
        scrollbarG.config(command=self.listbox.yview)
        button = Button(top, text="Calculate!", fg="red", font=("Courier", 22), command=self.ancestralNet)
        button.grid(row=1, column=0, columnspan=5, sticky=E+W)
        
    def ancestralNet(self): #Passess selected query and networks to ancestralNetwork.start function
        selected = []
        for i in self.listbox.curselection():
            selected.append(self.queryOrder.index(self.queriesAC[int(i)][1]))
        ancestralNetwork.start(selected)
            
######################################################################## 
####################################Settings menu items
                
####################Change database
    def UseDatabaseFile(self):
        database = askopenfilename(filetypes=[("Expression matrix",".exp")])
        if database !='':
            self.datafile = database.split('/')[len(database.split('/'))-1]
            configstring ="%s" % self.datafile
            DB = open(self.datafile, 'r').readlines()    
            self.LabVar.set("Using database: %s, Columns=%d, rows=%d." %(self.datafile.strip(), len(DB[1].split()), len(DB)))
            newconfig = open('dbconf.txt', 'w')
            newconfig.writelines(configstring)
            newconfig.close()
            self.printer("Using %s database\n" %(self.datafile))            
            try:
                config = open('dbconf.txt', 'r').readlines()
                self.currentDB = codecs.open('%s' % (config[0].rstrip()), mode='r', encoding='ASCII', errors='ignore').readlines()
                self.printer("Currently using database: %s consisting of %d probesets and %d chips.\n" % (config[0].rstrip(), len(self.currentDB), (len(self.currentDB[0].split())-3)))
                self.LabVar.set("Using database: %s, Columns=%d, rows=%d." %(config[0].rstrip(), len(self.currentDB[0].split()), len(self.currentDB)))
                nominator = []
                denominator = []
                self.annoDict = {}
                for i in range(len(self.currentDB)):
                    temp = []
                    splitta = self.currentDB[i].rstrip().split("\t")
                    self.annoDict[i] = splitta[0]+"\t"+splitta[1]+"\t"+splitta[2].replace("*"," ")+"\t"+splitta[3]+"\t"
                    for j in range(5,len(splitta)):
                        temp.append(float(splitta[j]))                    
                    expVector = numpy.array(temp)
                    nomi = expVector-(numpy.sum(expVector)/len(expVector))
                    denomi = numpy.sqrt(numpy.sum(nomi**2))
                    nominator.append(nomi)
                    denominator.append(denomi)
                self.nominator = numpy.array(nominator)
                self.denominator = numpy.array(denominator)
            except IOError:
                self.printer("Could not open current database.\n")
            
    def About(self):
        top = Toplevel()
        top.title("PlaNet standalone")
        label = Label(top, text="PlaNet standalone\nMarek Mutwil\ncontact: mr_mutwil@hotmail.com", width=30)
        OkButton = Button(top, text="OK", command=top.destroy)
        label.pack()
        OkButton.pack()
        
    def printer(self, text):   #prints text to message box
        #self.messenger.delete(1.0, END)
        self.messenger.insert(END, text)
    def clear(self):    #clears message box
        self.messenger.delete(1.0, END)
        
root = Tk()
root.title("Stand alone version of PlaNet")
app = App(root)
root.mainloop()
