#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dPSOgenerating as dPSO
import dataSaving as ds


#The model parameters must have type 'int' or 'float'. Elements of NumPy arrays ('numpy.int64' or 'numpy.float64') cannot be used.


# Generate one graph #####################################################################################################################################

directoryName_1ex='ExampleGraph' #all data will be saved in this directory
graphName_1ex = '_example' #this will be the graph's string identifier to be used in file names

#initialize the dPSOparameters class that will store the input parameters of the model, then set the parameters to the desired values:
iP = ds.dPSOparameters()
iP.d = 3 #generate networks in the d-dimensional hyperbolic space; 2<=d integer
iP.zeta = 1 #zeta=sqrt(-K), where K<0 is the curvature of the hyperbolic space; 0<zeta (usually zeta is set to 1)
iP.N = 100 #number of nodes; positive integer
iP.m = 2 #number of links created in each timestep after the mth step; positive integer; the average node degree will be approximately 2*m
iP.beta = 2/3 #popularity fading parameter; 0<beta<=1 (The absolute value of the exponent of the degree distribution is gamma=1+1/((d-1)*beta).)
iP.T = 0 #temperature; 0<=T<1/(d-1) or 1/(d-1)<T (T=0: the average clustering coefficient is maximized)
iP.numberOfGraphs = 1 #create one graph
iP.graphID = graphName_1ex #set the string identifier of the graph to be used in file names

#generate, save and return the graph, its parameters and the node coordinates:
[Rlist,GraphList,radialCoordList,CartCoordList] = dPSO.dPSO(inputParameters=iP,directoryName=directoryName_1ex,numORan='an',isWeighted=False) #the directory named directoryName_1ex will be created here, if it has not already been created
    #if numORan is 'num', then the cutoff distance R of the connection probability is calculated for each node by solving numerically Equation (S2.9) of the Supplementary Information, which expresses that the expected number of nodes connecting to the new node at its arrival is equal to m
    #if numORan is 'an', then the cutoff distance R of the connection probability is calculated faster, using the approximating formulas given by Equations (S2.11) and (S2.16) of the Supplementary Information
    #If isWeighted is True, then a weight calculated as 1/(1+hyperbolic distance between the two nodes in question) is assigned to each link at the end of the network growth. Otherwise, no weight is assigned to the links.
#the list with the cutoff distance of the connection probability in each time step is Rlist[0]
    #Rlist[0][j] is the cutoff distance at the appearance of the jth node (at time j+1)
    #Note that in the case of T=0 Rlist is just an empty list, since it is not needed for the network generation, as each node connects to the m number of hyperbolically closest nodes at T=0.
#the created NetworkX Graph is GraphList[0]
#the radial coordinates of all the network nodes in the native representation of the hyperbolic space are stored by the NumPy array radialCoordList[0]
    #radialCoordList[0][j] is the radial coordinate of the jth node of the network growth (j=0,1,...,N-1)
#the Cartesian coordinates of all the network nodes in the native representation of the hyperbolic space are stored by the N-by-d NumPy array CartCoordList[0]
    #CartCoordList[0][j,k] is the kth Cartesian coordinate of the jth node of the network growth (k=0,1,...,d-1 and j=0,1,...,N-1)



#load all the data from the txt files corresponding to the graph identifier graphName_1ex in the directory named directoryName_1ex:
[iP_loaded,Rlist_loaded,GraphList_loaded,radialCoordList_loaded,CartCoordList_loaded] = ds.load_dPSOgraph(directoryName=directoryName_1ex,isWeighted=False,whichGraph=graphName_1ex) #Note that in the case of T=0, Rlist_loaded is just an empty list, since it is not needed for the network generation, as each node connects to the m number of hyperbolically closest nodes at T=0.
#Here you can do what you want with the list Rlist_loaded[0], the NetworkX Graph GraphList_loaded[0], the NumPy array radialCoordList_loaded[0] and/or the N-by-d NumPy array CartCoordList_loaded[0].

##########################################################################################################################################################








#There are two options to generate a list of graphs using the SAME input parameters:

#1: Create (and save) one common dPSOparameters object (storing the model parameters) and set its numberOfGraphs parameter to the desired number of graphs. The graphs will be indexed by the 0, 1, ..., numberOfGraphs-1 integers. Save and load the graphs together, as the parts of a list. 

#2: Create (and save) a dPSOparameters object (containing the same model parameters) with a unique identifier (file name) for each graph and set the numberOfGraphs parameter to 1 in each case. Save and load the graphs one by one - in this way the overallocation of the memory resulting from storing all the graphs simultaneously in a list can be avoided. 



# Generate a list of graphs with the SAME input parameters - version 1 ###################################################################################

directoryName_exList1='ExampleGraphList1' #all data will be saved in this directory

#initialize the dPSOparameters class that will store the input parameters of the model, then set the parameters to the desired values:
iP = ds.dPSOparameters()
iP.d = 3 #generate networks in the d-dimensional hyperbolic space; 2<=d integer
iP.zeta = 1 #zeta=sqrt(-K), where K<0 is the curvature of the hyperbolic space; 0<zeta (usually zeta is set to 1)
iP.N = 100 #number of nodes; positive integer
iP.m = 2 #number of links created in each timestep after the mth step; positive integer; the average node degree will be approximately 2*m
iP.beta = 2/3 #popularity fading parameter; 0<beta<=1 (The absolute value of the exponent of the degree distribution is gamma=1+1/((d-1)*beta))
iP.T = 0 #temperature; 0<=T<1/(d-1) or 1/(d-1)<T (T=0: the average clustering coefficient is maximized)
iP.numberOfGraphs = 5 #positive integer
#The string identifier of the graph can be customized only if 1 graph is generated! For 1<iP.numberOfGraphs the graphs will be automatically indexed by the integers 0, 1, ..., iP.numberOfGraphs-1.

#generate, save and return the graphs, their parameters and the node coordinates:
[Rlist,GraphList,radialCoordList,CartCoordList] = dPSO.dPSO(inputParameters=iP,directoryName=directoryName_exList1,numORan='an',isWeighted=False) #the directory named directoryName_exList1 will be created here, if it has not already been created
    #if numORan is 'num', then the cutoff distance R of the connection probability is calculated for each node by solving numerically Equation (S2.9) of the Supplementary Information, which expresses that the expected number of nodes connecting to the new node at its arrival is equal to m
    #if numORan is 'an', then the cutoff distance R of the connection probability is calculated faster, using the approximating formulas given by Equations (S2.11) and (S2.16) of the Supplementary Information
    #If isWeighted is True, then a weight calculated as 1/(1+hyperbolic distance between the two nodes in question) is assigned to each link at the end of the network growth. Otherwise, no weight is assigned to the links.
#Rlist[g] is the list with the cutoff distance of the connection probability in each step of the growth of the gth graph
    #Rlist[g][j] is the cutoff distance at the appearance of the jth node (at time j+1) in the gth graph
    #Note that in the case of T=0 Rlist is just an empty list, since it is not needed for the network generation.
#GraphList[g] is the gth NetworkX Graph created by the d-dimensional E-PSO model (g=0,1,...,iP.numberOfGraphs-1)
#the radial coordinates of all the nodes of the gth graph are stored by the NumPy array radialCoordList[g]
    #radialCoordList[g][j] is the radial coordinate of the jth node of the network growth (j=0,1,...,N-1)
#the Cartesian coordinates of all the nodes of the gth graph are stored by the N-by-d NumPy array CartCoordList[g]
    #CartCoordList[g][j,k] is the kth Cartesian coordinate of the jth node of the network growth (k=0,1,...,d-1 and j=0,1,...,N-1)



#load all the data from the txt files in the directory named directoryName_exList1 with one function call:
[iP_loaded,Rlist_loaded,GraphList_loaded,radialCoordList_loaded,CartCoordList_loaded] = ds.load_dPSOgraph(directoryName=directoryName_exList1,isWeighted=False) #Note that in the case of T=0 Rlist_loaded is just an empty list, since it is not needed for the network generation.
#for g in range(iP_loaded.numberOfGraphs):
    #Here you can do what you want with the list Rlist_loaded[g], the NetworkX Graph GraphList_loaded[g], the NumPy array radialCoordList_loaded[g] and/or the N-by-d NumPy array CartCoordList_loaded[g].

##########################################################################################################################################################





# Generate a list of graphs with the SAME input parameters - version 2 ###################################################################################

directoryName_exList2='ExampleGraphList2' #all data will be saved in this directory
numOfGraphs = 5 #the number of graphs to be created
graphNameList = ["_"+str(i) for i in range(numOfGraphs)] #these will be the string identifiers of the graphs to be used in file names

#initialize the dPSOparameters class that will store the input parameters of the model, then set the parameters to the desired values:
iP = ds.dPSOparameters()
iP.d = 3 #generate networks in the d-dimensional hyperbolic space; 2<=d integer
iP.zeta = 1 #zeta=sqrt(-K), where K<0 is the curvature of the hyperbolic space; 0<zeta (usually zeta is set to 1)
iP.N = 100 #number of nodes; positive integer
iP.m = 2 #number of links created in each timestep after the mth step; positive integer; the average node degree will be approximately 2*m
iP.beta = 2/3 #popularity fading parameter; 0<beta<=1 (The absolute value of the exponent of the degree distribution is gamma=1+1/((d-1)*beta))
iP.T = 0 #temperature; 0<=T<1/(d-1) or 1/(d-1)<T (T=0: the average clustering coefficient is maximized)
iP.numberOfGraphs = 1 #create the graphs one by one

for g in range(numOfGraphs):
    iP.graphID = graphNameList[g] #set the string identifier of the given graph to be used in file names

    #generate and save the given graph and the node coordinates:
    [Rlist,GraphList,radialCoordList,CartCoordList] = dPSO.dPSO(inputParameters=iP,directoryName=directoryName_exList2,numORan='an',isWeighted=False) #the directory named directoryName_exList2 will be created here, if it has not already been created
        #if numORan is 'num', then the cutoff distance R of the connection probability is calculated for each node by solving numerically Equation (S2.9) of the Supplementary Information, which expresses that the expected number of nodes connecting to the new node at its arrival is equal to m
        #if numORan is 'an', then the cutoff distance R of the connection probability is calculated faster, using the approximating formulas given by Equations (S2.11) and (S2.16) of the Supplementary Information
        #If isWeighted is True, then a weight calculated as 1/(1+hyperbolic distance between the two nodes in question) is assigned to each link at the end of the network growth. Otherwise, no weight is assigned to the links.
    #the list with the cutoff distance of the connection probability in each time step is Rlist[0]
        #Rlist[0][j] is the cutoff distance at the appearance of the jth node (at time j+1)
        #Note that in the case of T=0 Rlist is just an empty list, since it is not needed for the network generation.
    #the created NetworkX Graph is GraphList[0]
    #the radial coordinates of all the network nodes in the native representation of the hyperbolic space is stored by the NumPy array radialCoordList[0]
        #radialCoordList[0][j] is the radial coordinate of the jth node of the network growth (j=0,1,...,N-1)
    #the Cartesian coordinates of all the network nodes in the native representation of the hyperbolic space is stored by the N-by-d NumPy array CartCoordList[0]
        #CartCoordList[0][j,k] is the kth Cartesian coordinate of the jth node of the network growth (k=0,1,...,d-1 and j=0,1,...,N-1)



#load all the data from the txt files in the directory named directoryName_exList2 one by one:
for g in range(numOfGraphs):
    [iP_loaded,Rlist_loaded,GraphList_loaded,radialCoordList_loaded,CartCoordList_loaded] = ds.load_dPSOgraph(directoryName=directoryName_exList2,isWeighted=False,whichGraph=graphNameList[g]) #Note that in the case of T=0 Rlist_loaded is just an empty list, since it is not needed for the network generation.
    #Here you can do what you want with the list Rlist_loaded[0], the NetworkX Graph GraphList_loaded[0], the NumPy array radialCoordList_loaded[0] and/or the N-by-d NumPy array CartCoordList_loaded[0].

