#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import networkx as nx
import datetime


#creating a directory for all of the files
def createDirectory(directoryName):
    filePath=os.getcwd()+"/"+directoryName+"/"
    directory=os.path.dirname(filePath)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)


#class for storing the input parameters of the 3-dimensional nPSO model together
class d3nPSOparameters:
    #create an object (named for example inputParam) for all of the input parameters with the command
        # inputParam=d3nPSOparameters()
    #this calling will set the parameters to the initial values given below
    #later the value of the parameters can be changed like inputParam.N=1000
    def __init__(self):
        self.zeta = 0 #zeta=sqrt(-K), where K<0 is the curvature of the hyperbolic space; 0<zeta (usually zeta is set to 1)
        self.N = 0 #number of nodes; positive integer
        self.m = 0 #number of links created in each timestep after the mth step; positive integer; the average node degree will be approximately 2*m
        self.beta = -1.0 #popularity fading parameter; 0<beta<=1 (The absolute value of the exponent of the degree distribution is gamma=1+1/(2*beta).)
        self.T = -1.0 #temperature; 0<=T<0.5 or 0.5<T (T=0: the average clustering coefficient is maximized.)
        self.numOfPeaks = 1 #number of peaks in the angular distribution of the nodes (i.e., the number of planted communities); positive integer
        self.probOfPeaks = [] #A list with with numOfPeaks number of elements, each of which is equal to the probability that the angular coordinates of a node will be sampled from the von Mises–Fisher distribution corresponding to a given peak; the sum of the elements must be 1
        self.muValues = [[]] #A list consisting of numOfPeaks number of lists of 3 elements, where each list of 3 elements corresponds to the Cartesian unit vector of the mean direction of a given peak in the angular distribution
        self.kappaValues = [] #a list with numOfPeaks number of real positive elements, each of which corresponding to the concentration of a given peak (the greater the value of kappa, the higher the concentration of the distribution around the mean direction mu; a smaller kappa leads to a broader angular patch; kappa=0 pertains to the uniform angular distribution of the nodes)
        self.numberOfGraphs = 0 #positive integer
        self.graphID = '' #the string identifier of the graph to be used in file names if numberOfGraphs==1
                          #if 1<numberOfGraphs: the graphs will be indexed by 0,1,...,numberOfGraphs-1


    #Function for SAVING the input parameters of the graphs into the directory named directoryName
    #the input parameters of a graph or graph list can be saved from a d3nPSOparameters object (named for example inputParam) with the command
        # inputParam.save(nameOfTheSavingDirectory)
    def save(self,directoryName):
        if self.numberOfGraphs==1:
            #save the input parameters in a readable way:
            ReportFileHandler_readable = open(os.getcwd() + "/" + directoryName + "/d3nPSOparameters_readable"+self.graphID+".txt", 'w')
            #create a file, which is more appropriate for the loading functions:
            ReportFileHandler = open(os.getcwd() + "/" + directoryName + "/d3nPSOparameters"+self.graphID+".txt", 'w')
        elif 1<self.numberOfGraphs:
            if self.graphID!='':
                print('Error: The graph identifier can be customized only if 1 graph is generated!')
            else:
                #save the input parameters in a readable way:
                ReportFileHandler_readable = open(os.getcwd() + "/" + directoryName + "/d3nPSOparameters_readable.txt", 'w')
                #create a file, which is more appropriate for the loading functions:
                ReportFileHandler = open(os.getcwd() + "/" + directoryName + "/d3nPSOparameters.txt", 'w')
        else:
            print('Error: numberOfGraphs can not be smaller than 1!')

        ReportFileHandler_readable.write('zeta:\t' + str(self.zeta))
        ReportFileHandler_readable.write('\nN:\t' + str(self.N))
        ReportFileHandler_readable.write('\nm:\t' + str(self.m))
        ReportFileHandler_readable.write('\nbeta:\t' + str(self.beta))
        ReportFileHandler_readable.write('\nT:\t' + str(self.T))
        ReportFileHandler_readable.write('\nnumber of graphs:\t' + str(self.numberOfGraphs))
        ReportFileHandler_readable.write('\nnumber of peaks:\t' + str(self.numOfPeaks))
        string_prob=str(self.probOfPeaks[0])
        for p in range(1,self.numOfPeaks):
            string_prob = string_prob+'\t'+str(self.probOfPeaks[p])
        ReportFileHandler_readable.write('\nprobability of choosing each peak:\t' + string_prob)
        ReportFileHandler_readable.write('\nCartesian coordinates of the mean direction vectors:\n')
        for p in range(self.numOfPeaks):
            ReportFileHandler_readable.write(str(self.muValues[p][0])+'\t'+str(self.muValues[p][1])+'\t'+str(self.muValues[p][2])+'\n')
        string_conc=str(self.kappaValues[0])
        for p in range(1,self.numOfPeaks):
            string_conc = string_conc+'\t'+str(self.kappaValues[p])
        ReportFileHandler_readable.write('concentration of each peak:\t' + string_conc)
        ReportFileHandler_readable.write('\n\n' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ReportFileHandler_readable.close()

        ReportFileHandler.write(str(self.zeta))
        ReportFileHandler.write("\n" + str(self.N))
        ReportFileHandler.write("\n" + str(self.m))
        ReportFileHandler.write("\n" + str(self.beta))
        ReportFileHandler.write("\n" + str(self.T))
        ReportFileHandler.write("\n" + str(self.numberOfGraphs))
        ReportFileHandler.write("\n" + str(self.numOfPeaks))
        string_prob=str(self.probOfPeaks[0])
        for p in range(1,self.numOfPeaks):
            string_prob = string_prob+'\t'+str(self.probOfPeaks[p])
        ReportFileHandler.write("\n" + string_prob)
        for p in range(self.numOfPeaks):
            ReportFileHandler.write("\n" + str(self.muValues[p][0])+'\t'+str(self.muValues[p][1])+'\t'+str(self.muValues[p][2]))
        string_conc='\n'+str(self.kappaValues[0])
        for p in range(1,self.numOfPeaks):
            string_conc = string_conc+'\t'+str(self.kappaValues[p])
        ReportFileHandler.write(string_conc)
        ReportFileHandler.close()


    #Function for LOADING the input parameters of graphs from the directory named directoryName
    #the input parameters of a previously saved graph (with a string identifier IDofTheGraph) or graph list can be loaded into a d3nPSOparameters object (named for example inputParam) with the command
        # inputParam.load(nameOfTheSavingDirectory,IDofTheGraph)
    def load(self,directoryName,graph_ID=''):
        if graph_ID != '': #load the parameters of one graph denoted by a previously given graph identifier
            self.graphID = graph_ID
            ReportFileHandler=open(os.getcwd()+"/"+directoryName+"/d3nPSOparameters"+self.graphID+".txt",'r')
            listOfLines=ReportFileHandler.readlines()
            ReportFileHandler.close()
            self.zeta=float(listOfLines[0])
            self.N=int(listOfLines[1])
            self.m=int(listOfLines[2])
            self.beta=float(listOfLines[3])
            self.T=float(listOfLines[4])
            self.numberOfGraphs=1

            self.numOfPeaks=int(listOfLines[6])

            self.probOfPeaks=[]
            listOfWords = listOfLines[7].split('\t')
            for j in range(self.numOfPeaks):
                self.probOfPeaks.append(float(listOfWords[j]))

            self.muValues=[]
            peakID = 0
            self.muValues.append([])
            for line in listOfLines[8:8+self.numOfPeaks]:
                listOfWords = line.split('\t')
                for j in range(3):
                    self.muValues[peakID].append(float(listOfWords[j]))
                peakID = peakID+1
                self.muValues.append([])

            self.kappaValues=[]
            listOfWords = listOfLines[8+self.numOfPeaks].split('\t')
            for j in range(self.numOfPeaks):
                self.kappaValues.append(float(listOfWords[j]))

        else: #default: self.graphID remains ''  ->  load the parameters of a graph list saved in the directory named directoryName
            ReportFileHandler=open(os.getcwd()+"/"+directoryName+"/d3nPSOparameters.txt",'r')
            listOfLines=ReportFileHandler.readlines()
            ReportFileHandler.close()
            self.zeta=float(listOfLines[0])
            self.N=int(listOfLines[1])
            self.m=int(listOfLines[2])
            self.beta=float(listOfLines[3])
            self.T=float(listOfLines[4])
            self.numberOfGraphs=int(listOfLines[5])

            self.numOfPeaks=int(listOfLines[6])

            self.probOfPeaks=[]
            listOfWords = listOfLines[7].split('\t')
            for j in range(self.numOfPeaks):
                self.probOfPeaks.append(float(listOfWords[j]))

            self.muValues=[]
            peakID = 0
            self.muValues.append([])
            for line in listOfLines[8:8+self.numOfPeaks]:
                listOfWords = line.split('\t')
                for j in range(3):
                    self.muValues[peakID].append(float(listOfWords[j]))
                peakID = peakID+1
                self.muValues.append([])

            self.kappaValues=[]
            listOfWords = listOfLines[8+self.numOfPeaks].split('\t')
            for j in range(self.numOfPeaks):
                self.kappaValues.append(float(listOfWords[j]))





####### Functions for data SAVING ########################################################################


#A function for saving the cutoff distance R of the connection probability in each time step. When more graphs are saved, one file is created for each graph. Each file contains only one colum of numbers.
#The file name always contains a graph identifier (if only one graph is saved, the identifier is the one given by the user in the corresponding d3nPSOparameters object; for a graph list each file name contains one of the integers 0,1,2,...,numberOfGraphs-1).

#inputParameters=d3nPSOparameters object with the input parameters of the model
#listToSave=list of lists, where each list contains the cutoff distance of the connection probability in each time step for a given graph. (In one listToSave object all the corresponding graphs' input parameters are given by the same d3nPSOparameters object.) When saving only one list named givenList, set the argument listToSave to [givenList].
#fileName is the beginning of the file names. It can be for example "cutoffDistances" for saving the cutoff distance of the connection probability in each time step.
#directoryName=the name of the directory in which the file(s) will be saved
def save_Rlist(inputParameters,listToSave,fileName,directoryName):
    if inputParameters.numberOfGraphs==1: #listToSave is a list containing only one list
        wholeFileName=fileName+inputParameters.graphID+".txt"
        #print("Save a list to "+wholeFileName)
        fileHandler=open(os.getcwd()+"/"+directoryName+"/"+wholeFileName,"w")
        for j in range(inputParameters.N):
            fileHandler.write(str(listToSave[0][j]) + "\n") #list element for the jth node
        fileHandler.close()

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #listToSave is a list containing more than one lists
        for i in range(inputParameters.numberOfGraphs):
            wholeFileName=fileName+"_"+str(i)+".txt"
            #print("Save a list to "+wholeFileName)
            fileHandler=open(os.getcwd()+"/"+directoryName+"/"+wholeFileName,"w")
            for j in range(inputParameters.N):
                fileHandler.write(str(listToSave[i][j])+"\n") #list element for the jth node
            fileHandler.close()
    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only if 1 graph is generated!')




#Save the edge list of the created graph(s) in txt file(s). When more graphs are saved, one file is created for each graph.
#The file name always contains a graph identifier (if one graph is saved, the identifier is the one given by the user in the corresponding d3nPSOparameters object; for a graph list each file name contains one of the integers 0,1,2,...,numberOfGraphs-1).
#If isWeighted==False, than each line is of the form
#id1	id2
#if node id1 is connected with node id2.
#If isWeighted==True, than each line is of the form
#id1	id2	weight
#if node id1 is connected with node id2 with an edge.
#The node identifiers are integers (0,1,2,...,N-1). The values within each row are separated by tabs.

#inputParameters=d3nPSOparameters object with the input parameters of the model
#GraphList=A list of the created NetworkX Graphs. (In one list all the graphs' input parameters are given by the same d3nPSOparameters object.) When saving only one graph (e.g. G), set the argument GraphList to [G].
#isWeighted=False, if the edges are not weighted and True otherwise
#directoryName=the name of the directory in which the file(s) will be saved
def saveGraphs_EdgeList(inputParameters,GraphList,isWeighted,directoryName):
    if inputParameters.numberOfGraphs==1: #GraphList is a list containing only one NetworkX Graph
        edgeListFileName="edgeList"+inputParameters.graphID+".txt"
        #print("Save graph "+inputParameters.graphID+"'s edge list in txt: "+edgeListFileName)
        if isWeighted:
            nx.write_edgelist(GraphList[0],os.getcwd()+"/"+directoryName+"/"+edgeListFileName,delimiter='\t',data=['weight'])
        else: #not weighted edges
            nx.write_edgelist(GraphList[0],os.getcwd()+"/"+directoryName+"/"+edgeListFileName,delimiter='\t',data=False)

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #GraphList is a list containing more than one NetworkX Graphs
        if isWeighted:
            for i in range(inputParameters.numberOfGraphs):
                edgeListFileName="edgeList_"+str(i)+".txt"
                #print("Save graph "+str(i)+"'s edge list in txt: "+edgeListFileName)
                nx.write_edgelist(GraphList[i],os.getcwd()+"/"+directoryName+"/"+edgeListFileName,delimiter='\t',data=['weight'])
        else: #not weighted edges
            for i in range(inputParameters.numberOfGraphs):
                edgeListFileName="edgeList_"+str(i)+".txt"
                #print("Save graph "+str(i)+"'s edge list in txt: "+edgeListFileName)
                nx.write_edgelist(GraphList[i],os.getcwd()+"/"+directoryName+"/"+edgeListFileName,delimiter='\t',data=False)
    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only if 1 graph is generated!')




#Save the radial coordinates and the Cartesian coordinates of the network nodes in the native representation of the hyperbolic space. When more graphs are saved, one file is created for each graph.
#The file name always contains a graph identifier (if only one graph is saved, the identifier is the one given by the user in the corresponding d3nPSOparameters object; for a graph list each file name contains one of the integers 0,1,2,...,numberOfGraphs-1).
#Each line is of the form
#id	radialCoord	CartCoord1	CartCoord2	...
#where id is the node identifier (0<=id<=numberOfNodes-1 integer), radialCoord is the radial coordinate of this node and CartCoord1, CartCoord2, etc. are the Cartesian coordinates of this node. The values within each row are separated by tabs.

#inputParameters=d3nPSOparameters object with the input parameters of the model
#radialCoordList,CartCoordList: Both are lists of NumPy arrays with the radial and Cartesian node coordinates. (In one list all the corresponding graphs' input parameters are given by the same d3nPSOparameters object.) When saving only one array for each coordinate type (e.g. r and Coords), set the arguments radialCoordList to [r] and CartCoordList to [Coords].
#fileName is the beginning of the file names. It can be for example
    #"originalCoordinates" for saving the original coordinates of the nodes generated by the d-dimensional PSO model
    #"embeddedCoordinates" for saving the coordinates produced by an embedding method
#directoryName=the name of the directory in which the file(s) will be saved
def saveCoordinates(inputParameters,radialCoordList,CartCoordList,fileName,directoryName):
    if inputParameters.numberOfGraphs==1: #radialCoordList is a list containing only one NumPy array and CartCoordList is a list containing only one NumPy array
        coordFileName=fileName+inputParameters.graphID+".txt"
        #print("Save coordinate arrays to "+coordFileName)
        fileHandler=open(os.getcwd()+"/"+directoryName+"/"+coordFileName,"w")
        for j in range(inputParameters.N):
            fileHandler.write(str(j) + "\t") #node identifier
            fileHandler.write(str(radialCoordList[0][j])+"\t") #radial coordinate of the jth node
            for k in range(3-1):
                fileHandler.write(str(CartCoordList[0][j,k])+"\t") #kth Cartesian coordinate of the jth node
            fileHandler.write(str(CartCoordList[0][j,3-1])+"\n") #the last Cartesian coordinate of the jth node
        fileHandler.close()

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #radialCoordList is a list containing more than one NumPy arrays and CartCoordList is a list containing more than one NumPy arrays
        for i in range(inputParameters.numberOfGraphs):
            coordFileName=fileName+"_"+str(i)+".txt"
            #print("Save coordinate arrays to "+coordFileName)
            fileHandler=open(os.getcwd()+"/"+directoryName+"/"+coordFileName,"w")
            for j in range(inputParameters.N):
                fileHandler.write(str(j) + "\t") #node identifier
                fileHandler.write(str(radialCoordList[i][j])+"\t") #radial coordinate of the jth node
                for k in range(3-1):
                    fileHandler.write(str(CartCoordList[i][j,k])+"\t") #kth Cartesian coordinate of the jth node
                fileHandler.write(str(CartCoordList[i][j,3-1])+"\n") #the last Cartesian coordinate of the jth node
            fileHandler.close()
    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only if 1 graph is generated!')



def savePlantedCommunities(inputParameters,angularGroupIDlist,fileName,directoryName):
    if inputParameters.numberOfGraphs==1: #angularGroupIDlist is a list containing only one list of community indexes
        commStructFileName=fileName+inputParameters.graphID+".txt"
        fileHandler=open(os.getcwd()+"/"+directoryName+"/"+commStructFileName,"w")
        for j in range(inputParameters.N):
            fileHandler.write(str(j) + "\t") #node identifier
            fileHandler.write(str(angularGroupIDlist[0][j])+"\n") #group index
        fileHandler.close()

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #angularGroupIDlist is a list containing more than one list of community indexes
        for i in range(inputParameters.numberOfGraphs):
            commStructFileName=fileName+"_"+str(i)+".txt"
            fileHandler=open(os.getcwd()+"/"+directoryName+"/"+commStructFileName,"w")
            for j in range(inputParameters.N):
                fileHandler.write(str(j) + "\t") #node identifier
                fileHandler.write(str(angularGroupIDlist[i][j])+"\n") #group index
            fileHandler.close()
    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only if 1 graph is generated!')



#Save all data regarding the generated d-dimensional PSO network(s) in the already existing directory named directoryName

#inputParameters=d3nPSOparameters object with the input parameters of the model
#Rlist=a list of the lists containing the cutoff distance of the connection probability in each step of the growth of a graph (Rlist[i][j] is the cutoff distance at the appearance of node j in graph i); if inputParameters.T==0, no such list is created during the network growth, therefore this argument will be ignored here
#GraphList=a list of the created NetworkX Graphs
#isWeighted=False, if the edges are not weighted and True otherwise
#radialCoordList,CartCoordList=lists of NumPy arrays with the node coordinates (radial and Cartesian) in the native representation of the hyperbolic space
#angularGroupIDlist=list of NumPy arrays containing the identifier of the chosen peak of the angular distribution for each node of the given graph
#directoryName=the name of the directory in which the file(s) will be saved
def saveAll(inputParameters,Rlist,GraphList,isWeighted,radialCoordList,CartCoordList,angularGroupIDlist,directoryName):
    createDirectory(directoryName) #create the directory for storing all results
    inputParameters.save(directoryName)
    if inputParameters.T<0 or inputParameters.T==1/(3-1): #d=3
        print("Choose 0<=T<1/(d-1) or 1/(d-1)<T!")
    elif 0<inputParameters.T: #for T=0 the cutoff distance of the connection probability has no role in the network generation
        save_Rlist(inputParameters,Rlist,"cutoffDistances",directoryName)
    saveGraphs_EdgeList(inputParameters,GraphList,isWeighted,directoryName)
    saveCoordinates(inputParameters,radialCoordList,CartCoordList,"originalCoordinates",directoryName)
    savePlantedCommunities(inputParameters,angularGroupIDlist,"communities_angular",directoryName)






###### Functions for data LOADING ########################################################################


#Load the cutoff distance R of the connection probability in each time step from the corresponding txt file(s) in the directory named directoryName. Return a LIST of lists, even for only one graph.
#Each file must contain only one colum of numbers.

#directoryName=the name of the directory from which the file(s) will be loaded
#fileName is the beginning of the file name(s). It can be for example "cutoffDistances" for loading the cutoff distance of the connection probability in each time step.
#inputParameters=d3nPSOparameters object with the input parameters of the model
def load_Rlist(directoryName,fileName,inputParameters):
    if inputParameters.numberOfGraphs==1: #load the data for a given graph
        wholeFileName = fileName+inputParameters.graphID+".txt"
        #print("Load a list from: "+directoryName+"/"+wholeFileName)
        Rvalues=list(np.loadtxt(os.getcwd()+"/"+directoryName+"/"+wholeFileName))
        return [Rvalues] #return the list as the only element of a list

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #load the data for a list of NetworkX Graphs created using the parameters given by the same d3nPSOparameters object, i.e. by inputParameters
        Rlist = [] #a list containing the desired list for all graphs
        for i in range(inputParameters.numberOfGraphs):
            wholeFileName = fileName+"_"+str(i)+".txt"
            #print("Load a list from: "+directoryName+"/"+wholeFileName)
            Rvalues=list(np.loadtxt(os.getcwd()+"/"+directoryName+"/"+wholeFileName))
            Rlist.append(Rvalues)
        return Rlist #return the list of lists

    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only for 1 graph!')




#Load d-dimensional PSO network(s) from txt file(s) containing its/their edge list(s) in the directory named directoryName. Return a LIST of NetworkX Graphs, even if only one graph is loaded.
#Each line must be of the form
#id1	id2
#or
#id1	id2	weight
#if node id1 is connected with node id2 with an edge.
#The node identifiers must be integers (0,1,2,...,N-1). The values within each row must be separated by tabs.

#directoryName=the name of the directory from which the file(s) will be loaded
#inputParameters=d3nPSOparameters object with the input parameters of the model
#if isWeighted==True, than weights will be assigned to the edges, otherwise not
def loadGraphList_edgeList(directoryName,inputParameters,isWeighted):
    if inputParameters.numberOfGraphs==1: #load a given graph
        edgeListFileName="edgeList"+inputParameters.graphID+".txt"
        #print("Load graph "+inputParameters.graphID+" from: "+directoryName+"/"+edgeListFileName)
        G = nx.Graph() #the initial empty network
        edgeList = [] #list with the edge tuples for the add_edges_from function
        fileHandler = open(os.getcwd()+"/"+directoryName+"/"+edgeListFileName,"r")
        listOfLines = fileHandler.readlines()
        fileHandler.close()
        if isWeighted: #load the edge weights
            for line in listOfLines:
                listOfWords = line.split('\t')
                node1 = int(listOfWords[0])
                node2 = int(listOfWords[1])
                weight = float(listOfWords[2])
                edgeList.append((node1,node2,weight))
            G.add_nodes_from(list(range(inputParameters.N))) #The node identifiers are always the integers 0, 1, 2, ..., N-1
            G.add_weighted_edges_from(edgeList)
            return [G] #return the desired NetworkX Graph as the only element of a list
        else: #discard the edge weights / the edge list consists of only 2 columns
            for line in listOfLines:
                listOfWords = line.split('\t')
                node1 = int(listOfWords[0])
                node2 = int(listOfWords[1])
                edgeList.append((node1, node2))
            G.add_nodes_from(list(range(inputParameters.N))) #The node identifiers are always the integers 0, 1, 2, ..., N-1
            G.add_edges_from(edgeList)
            return [G] #return the desired NetworkX Graph as the only element of a list

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #load a list of NetworkX Graphs created using the parameters given by the same d3nPSOparameters object, i.e. by inputParameters
        GraphList = [] #list of the graphs
        if isWeighted: #load the edge weights
            for i in range(inputParameters.numberOfGraphs):
                edgeListFileName=("edgeList_"+str(i)+".txt")
                #print("Load graph "+str(i)+" from: "+directoryName+"/"+edgeListFileName)
                G=nx.Graph() #the initial empty network
                edgeList=[] #list with the edge tuples for the add_edges_from function
                fileHandler=open(os.getcwd()+"/"+directoryName+"/"+edgeListFileName, "r")
                listOfLines = fileHandler.readlines()
                fileHandler.close()
                for line in listOfLines:
                    listOfWords=line.split('\t')
                    node1=int(listOfWords[0])
                    node2=int(listOfWords[1])
                    weight=float(listOfWords[2])
                    edgeList.append((node1,node2,weight))
                G.add_nodes_from(list(range(inputParameters.N))) #The node identifiers are always the integers 0, 1, 2, ..., N-1
                G.add_weighted_edges_from(edgeList)
                GraphList.append(G)
            return GraphList #return a list of the NetworkX Graphs
        else: #discard the edge weights / the edge list consists of only 2 columns
            for i in range(inputParameters.numberOfGraphs):
                edgeListFileName=("edgeList_"+str(i)+".txt")
                #print("Load graph "+str(i)+" from: "+directoryName+"/"+edgeListFileName)
                G=nx.Graph() #the initial empty network
                edgeList=[] #list with the edge tuples for the add_edges_from function
                fileHandler=open(os.getcwd()+"/"+directoryName+"/"+edgeListFileName, "r")
                listOfLines = fileHandler.readlines()
                fileHandler.close()
                for line in listOfLines:
                    listOfWords=line.split('\t')
                    node1=int(listOfWords[0])
                    node2=int(listOfWords[1])
                    edgeList.append((node1,node2))
                G.add_nodes_from(list(range(inputParameters.N))) #The node identifiers are always the integers 0, 1, 2, ..., N-1
                G.add_edges_from(edgeList)
                GraphList.append(G)
            return GraphList #return a list of the NetworkX Graphs

    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only for 1 graph!')




#Load the node coordinates in the native representation of the d-dimensional hyperbolic space from the corresponding txt file(s) in the directory named directoryName. Return a LIST of NumPy arrays containing the node coordinates, even if only one graph's coordinate set is loaded.
#Each line must be of the form
#id	radialCoord	CartCoord1	CartCoord2	...
#where id is the node identifier (0<=id<=numberOfNodes-1 integer), radialCoord is the radial coordinate of this node and CartCoord1, CartCoord2, etc. are the Cartesian coordinates of this node. The values within each row must be separated by tabs.

#directoryName=the name of the directory from which the file(s) will be loaded
#fileName is the beginning of the file name(s). It can be for example
    #"originalCoordinates" for loading the original coordinates of the nodes generated by the d-dimensional PSO model
    #"embeddedCoordinates" for loading the coordinates produced by an embedding method
#inputParameters=d3nPSOparameters object with the input parameters of the model
def loadCoordinates(directoryName,fileName,inputParameters):
    r=np.zeros(inputParameters.N) #radial coordinates of the nodes in one graph
    Coords=np.zeros((inputParameters.N,3)) #Cartesian coordinates of the nodes in one graph

    if inputParameters.numberOfGraphs==1: #load the data for a given graph
        coordFileName = fileName+inputParameters.graphID+".txt"
        #print("Load coordinate arrays from: "+directoryName+"/"+coordFileName)
        fileHandler = open(os.getcwd()+"/"+directoryName+"/"+coordFileName,"r")
        listOfLines = fileHandler.readlines()
        fileHandler.close()
        for line in listOfLines:
            listOfWords = line.split('\t')
            nodeID = int(listOfWords[0])
            r[nodeID] = float(listOfWords[1])
            for j in range(2,2+3):
                Coords[nodeID,j-2] = float(listOfWords[j])
        return [[r],[Coords]] #return the NumPy array in both cases as the only element of a list

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #load the data for a list of NetworkX Graphs created using the parameters given by the same d3nPSOparameters object, i.e. by inputParameters
        radialCoordList = [] #list of the radial coordinates of all graphs
        CartCoordList = [] #list of the Cartesian coordinates of all graphs
        for i in range(inputParameters.numberOfGraphs):
            coordFileName = fileName+"_"+str(i)+".txt"
            #print("Load coordinate arrays from: "+directoryName+"/"+coordFileName)
            fileHandler = open(os.getcwd()+"/"+directoryName+"/"+coordFileName,"r")
            listOfLines = fileHandler.readlines()
            fileHandler.close()
            for line in listOfLines:
                listOfWords = line.split('\t')
                nodeID = int(listOfWords[0])
                r[nodeID] = float(listOfWords[1])
                for j in range(2,2+3):
                    Coords[nodeID,j-2] = float(listOfWords[j])
            radialCoordList.append(r)
            CartCoordList.append(Coords)
        return [radialCoordList,CartCoordList] #return two lists of the NumPy arrays

    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only for 1 graph!')



def load_plantedCommunities(directoryName,fileName,inputParameters):
    angularGroupIDs = np.zeros(inputParameters.N)
    if inputParameters.numberOfGraphs==1: #load the data for a given graph
        wholeFileName = fileName+inputParameters.graphID+".txt"
        fileHandler = open(os.getcwd()+"/"+directoryName+"/"+wholeFileName,"r")
        listOfLines = fileHandler.readlines()
        fileHandler.close()
        for line in listOfLines:
            listOfWords = line.split('\t')
            nodeID = int(listOfWords[0])
            angularGroupIDs[nodeID] = int(listOfWords[1])
        return [list(angularGroupIDs)] #return the list as the only element of a list

    elif 1<inputParameters.numberOfGraphs and inputParameters.graphID=='': #load the data for a list of NetworkX Graphs created using the parameters given by the same d3nPSOparameters object, i.e. by inputParameters
        angularGroupIDlist = [] #a list containing the desired list for all graphs
        for i in range(inputParameters.numberOfGraphs):
            wholeFileName = fileName+"_"+str(i)+".txt"
            fileHandler = open(os.getcwd()+"/"+directoryName+"/"+wholeFileName,"r")
            listOfLines = fileHandler.readlines()
            fileHandler.close()
            for line in listOfLines:
                listOfWords = line.split('\t')
                nodeID = int(listOfWords[0])
                angularGroupIDs[nodeID] = int(listOfWords[1])
            angularGroupIDlist.append(angularGroupIDs)
        return angularGroupIDlist #return the list of lists

    else: #1<inputParameters.numberOfGraphs and inputParameters.graphID!=''
        print('Error: The graph identifier can be customized only for 1 graph!')



#Load all data regarding the previously generated 3-dimensional nPSO network(s) from the directory named directoryName

#if isWeighted==True, than weights will be assigned to the edges, otherwise not
#Use whichGraph='' (default) to load a whole graph list from the directory named directoryName. (In one list all the corresponding graphs' input parameters are given by the same d3nPSOparameters object.)
#To load a graph with a customized identifier (and thus with a customized d3nPSOparameters object), set whichGraph to the string identifier of the desired graph.
def load_d3nPSOgraph(directoryName,isWeighted,whichGraph=''):
    params = d3nPSOparameters() #create a d3nPSOparameters object for storing the model parameters
    params.load(directoryName,whichGraph) #load the input parameters of the d-dimensional PSO model, with which the network(s) was/were generated
    if params.T==0:
        Rlist = [] #in the case of T=0 Rlist is not needed for the network generation
    else: #0<params.T<1/(params.d-1) or 1/(params.d-1)<params.T
        Rlist = load_Rlist(directoryName,"cutoffDistances",params)
    Rlist = []
    GraphList = loadGraphList_edgeList(directoryName,params,isWeighted) #load a list of NetworkX Graphs
    [radialCoordList,CartCoordList]=loadCoordinates(directoryName,"originalCoordinates",params) #load two lists of NumPy arrays
    angularGroupIDlist=load_plantedCommunities(directoryName,"communities_angular",params) #load the planted community structure
    return [params, Rlist, GraphList, radialCoordList, CartCoordList,angularGroupIDlist]

