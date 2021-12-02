#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import d3nPSOdataSaving as ds
import numpy as np
import networkx as nx
import random as rand
import math
import matlab.engine
import os
from tensorflow_probability import distributions as dist #von Mises–Fisher distribution, i.e. isotropic Gaussian distribution defined over a 2-dimensional sphere (the highest allowed dimension of the sphere is 4)


#Generate random networks in the hyperbolic space with thed-dimensional PSO model and save the results (edgelist, node coordinates) in the directory named directoryName.
#inputParameters=dPSOparameters object (defined in dataSaving.py) with the input parameters of the model
#numORan is a string specifying what type of calculations will be used for determining the cutoff distance of the connection probability:
    #if numORan=='num', the cutoff distances will be calculated by solving the equation m=i*P(i) numerically
    #if numORan=='an', the cutoff distances will be calculated analytically, using approximating formulas
    #Note that the parameter numORan will not be used if T=0, since in this case the cutoff distance is not needed for the network generation.
#if isWeighted==True, the edges will be weighted according to the hyperbolic distances (w_ij=1/(1+h_ij))
#This functon returns a list of the generated NetworkX Graphs (GraphList), a list of the NumPy arrays storing the radial coordinates of the nodes in the native representation (radialCoordList) and a list of the N-by-d NumPy arrays storing all the Cartesian coordinates of the nodes in the native representation.
def d3nPSO(inputParameters,directoryName,numORan,isWeighted):
    probNorm = sum(inputParameters.probOfPeaks)
    if probNorm!=1:
        for i in range(inputParameters.numOfPeaks):
            inputParameters.probOfPeaks[i] = inputParameters.probOfPeaks[i]/probNorm #normalization
    print("probability of peaks after normalisation: "+str(inputParameters.probOfPeaks))
    for i in range(inputParameters.numOfPeaks):
        inputParameters.muValues[i] = list(np.array(inputParameters.muValues[i])/np.linalg.norm(inputParameters.muValues[i]))
    print("mean directions as unit vectors: "+str(inputParameters.muValues))
    if inputParameters.T==0: #deterministic edge formation
        [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist] = d3nPSO_T0(inputParameters.zeta, inputParameters.N, inputParameters.m, inputParameters.beta, inputParameters.probOfPeaks, inputParameters.muValues, inputParameters.kappaValues, inputParameters.numberOfGraphs, isWeighted) #note that Rlist is empty for T==0
        ds.saveAll(inputParameters,Rlist,GraphList,isWeighted,radialCoordList,CartCoordList,angularGroupIDlist,directoryName) #save the edge list and the node coordinates for all the generated networks
        return [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist]

    elif inputParameters.T<(1/(3-1)): #d=3
        [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist] = d3nPSO_Tnot0_belowLimit(inputParameters.zeta, inputParameters.N, inputParameters.m, inputParameters.beta, inputParameters.T, inputParameters.probOfPeaks, inputParameters.muValues, inputParameters.kappaValues, inputParameters.numberOfGraphs, numORan, isWeighted)
        ds.saveAll(inputParameters,Rlist,GraphList,isWeighted,radialCoordList,CartCoordList,angularGroupIDlist,directoryName) #save the edge list and the node coordinates for all the generated networks
        return [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist]

    elif (1/(3-1))<inputParameters.T: #d=3
        [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist] = d3nPSO_TaboveLimit(inputParameters.zeta, inputParameters.N, inputParameters.m, inputParameters.beta, inputParameters.T, inputParameters.probOfPeaks, inputParameters.muValues, inputParameters.kappaValues, inputParameters.numberOfGraphs, numORan, isWeighted)
        ds.saveAll(inputParameters,Rlist,GraphList,isWeighted,radialCoordList,CartCoordList,angularGroupIDlist,directoryName) #save the edge list and the node coordinates for all the generated networks
        return [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist]

    else:
        print("Choose 0<=T<0.5 or 0.5<T!")



def d3nPSO_T0(zeta,N,m,beta,probOfPeaks,muValues,kappaValues,numberOfGraphs,isWeighted):
    GraphList = [] #list of the generated NetworkX Graphs
    radialCoordList = [] #radialCoordList[g]=NumPy array with the radial coordinates of the gth graph's nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
    CartCoordList = [] #CartCoordList[g]=N-by-3 NumPy array containing the Cartesian coordinates describing the positions of the gth graph's nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
    angularGroupIDlist = []

    for g in range(numberOfGraphs):
        #Initialization
        G = nx.Graph() #the initial empty network
        r = np.zeros(N) #radial coordinates of the N nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
        r_initial = (2/zeta)*np.log(np.linspace(1,N,N)) #initial radial coordinates of the N nodes ((2/zeta)*ln(timestep)) in the native representation
        Coords = np.zeros((N,3)) #Cartesian coordinates of the N nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
        angularGroupIDs = -np.ones(N,dtype=int)

        #In each timestep one new node is added to the graph.
        #The time is i+1 when node i appears!!!
        #T=0: deterministic edge formation, i.e. the new node i connects to the m hyperbolically closest nodes

        #Add the first two nodes (indexed by 0 and 1) to the graph
        G.add_node(0)
        r[0] = r_initial[0] #set the current radial coordinate of node 0 to its initial value
        Coords[0,:] = np.zeros(3) #set the initial Cartesian coordinates of node 0
        G.add_node(1)
        r[1] = r_initial[1] #set the current radial coordinate of node 1 to its initial value
        [angularGroupIDs[1],Coords[1,:]] = generateCartCoords(r[1],probOfPeaks,muValues,kappaValues) #generate the initial Cartesian coordinates of node 1
        #simulate popularity fading, if needed
        if beta<1: #there is popularity fading
            r[0] = beta*r_initial[0]+(1-beta)*r_initial[1]
            [angularGroupIDs[0],Coords[0,:]] = generateCartCoords(r[0],probOfPeaks,muValues,kappaValues) #update the Cartesian coordinates of node 0
            #Note that initially the length of node 0's coordinate array is 0 (r[0]=0), therefore the below applied update rule Coords[0,:] = r_0_new*Coords[0,:]/r[0] would not work here.
        G.add_edges_from([(0,1)]) #add the new edge to the graph

        #Continue the network growth
        for i in range(2,N):
            print("Node "+str(i))
            G.add_node(i) #add the new node to the graph
            r[i] = r_initial[i] #set the current radial coordinate of node i to its initial value
            [angularGroupIDs[i],Coords[i,:]] = generateCartCoords(r[i],probOfPeaks,muValues,kappaValues) #generate the initial Cartesian coordinates of node i

            edgeList = [] #list of the edges (node identifier tuples) created at the appearance of node i
            #create m edges connecting node i to the previously appeared nodes
            #note that at time i+1 the number of previously appeared nodes is i, therefore the maximum number of new edges is i
            if i<=m: #connect the new node to all of the previously appeared nodes
                for j in range(i): #j=indices of the previously appeared nodes
                    edgeList.append((j,i))
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
            else: #m<i: connect the new node to the m hyperbolically closest previously appeared nodes
                hypDistances = np.zeros(i) #for the hyperbolic distances between the new node (i) and the already existing nodes
                for j in range(i): #j=indices of the existing nodes
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
                    hypDistances[j] = hypDist(zeta,Coords[i,:],Coords[j,:],r[i],r[j])
                closestNodeIndices = np.argpartition(hypDistances,m)[:m] #the indices of the m closest nodes
                for k in range(m):
                    edgeList.append((closestNodeIndices[k],i))
            G.add_edges_from(edgeList) #add the new edges to the graph

        #one graph has been created
        if isWeighted:
            weighting(zeta,G,Coords,r)
        print("Graph "+str(g)+" is ready.")
        GraphList.append(G)
        radialCoordList.append(r)
        CartCoordList.append(Coords)
        angularGroupIDlist.append(angularGroupIDs)
    #all graphs are ready
    return [[],GraphList,radialCoordList,CartCoordList,angularGroupIDlist] #[]: in the case of T=0 Rlist is not calculated, since it is not needed for the network generation



def d3nPSO_Tnot0_belowLimit(zeta,N,m,beta,T,probOfPeaks,muValues,kappaValues,numberOfGraphs,numORan,isWeighted):
    Rlist = [] #Rlist[g]=list with the cutoff distance of the connection probability in each step of the growth of the gth graph
    GraphList = [] #list of the generated NetworkX Graphs
    radialCoordList = [] #radialCoordList[g]=NumPy array with the radial coordinates of the gth graph's nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
    CartCoordList = [] #CartCoordList[g]=N-by-3 NumPy array containing the Cartesian coordinates describing the positions of the gth graph's nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
    angularGroupIDlist = []

    #Calculate the cutoff distance of the connection probability in each step of the network growth (the cutoff distance is the same for each graph)
    Rvalues = cutoffDists(zeta,N,m,beta,T,numORan) #the cutoff distance of the connection probability for node i is Rvalues[i]

    for g in range(numberOfGraphs):
        #Initialization
        G = nx.Graph() #the initial empty network
        r = np.zeros(N) #radial coordinates of the N nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
        r_initial = (2/zeta)*np.log(np.linspace(1,N,N)) #initial radial coordinates of the N nodes ((2/zeta)*ln(timestep)) in the native representation
        Coords = np.zeros((N,3)) #Cartesian coordinates of the N nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
        angularGroupIDs = -np.ones(N,dtype=int)

        #In each timestep one new node is added to the graph.
        #The time is i+1 when node i appears!!!
        #0<T: the new node i chooses its neighbors according to the probability distribution p_i[j]

        #Add the first two nodes (indexed by 0 and 1) to the graph
        G.add_node(0)
        r[0] = r_initial[0] #set the current radial coordinate of node 0 to its initial value
        Coords[0,:] = np.zeros(3) #set the initial Cartesian coordinates of node 0
        print("Node 0 does not create connections - the value of R_0 is irrelevant.")
        G.add_node(1)
        r[1] = r_initial[1] #set the current radial coordinate of node 1 to its initial value
        [angularGroupIDs[1],Coords[1,:]] = generateCartCoords(r[1],probOfPeaks,muValues,kappaValues) #generate the initial Cartesian coordinates of node 1
        #simulate popularity fading, if needed
        if beta<1: #there is popularity fading
            r[0] = beta*r_initial[0]+(1-beta)*r_initial[1]
            [angularGroupIDs[0],Coords[0,:]] = generateCartCoords(r[0],probOfPeaks,muValues,kappaValues) #update the Cartesian coordinates of node 0
            #Note that initially the length of node 0's coordinate array is 0 (r[0]=0), therefore the below applied update rule Coords[0,:] = r_0_new*Coords[0,:]/r[0] would not work here.
        #edge creation
        G.add_edges_from([(0,1)]) #add the new edge to the graph
        print("Connect node 1 to all of the previously appeared nodes - the value of R_1 is irrelevant.")

        #Continue the network growth
        for i in range(2,N):
            print("Node "+str(i))
            G.add_node(i) #add the new node to the graph
            r[i] = r_initial[i] #set the current radial coordinate of node i to its initial value
            [angularGroupIDs[i],Coords[i,:]] = generateCartCoords(r[i],probOfPeaks,muValues,kappaValues) #generate the Cartesian coordinates of node i

            edgeList = [] #list of the edges (node identifier tuples) created at the appearance of node i
            #create m edges connecting node i to the previously appeared nodes
            #note that at time i+1 the number of previously appeared nodes is i, therefore the maximum number of new edges is i
            if i<=m: #connect the new node to all of the previously appeared nodes
                print("Connect node "+str(i)+" to all of the previously appeared nodes - the value of R_"+str(i)+" is irrelevant.")
                for j in range(i): #j=indices of the previously appeared nodes
                    edgeList.append((j,i))
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
            else: #m<i: connect the new node to m previously appeared nodes chosen according to the probability distribution p_i[j]
                R_i = Rvalues[i] #set the cutoff distance of the connection probability
                #determine the connection probabilities
                hypDistances = np.zeros(i) #for the hyperbolic distances between the new node and the already existing nodes
                p_i = np.zeros(i) #for the connection probabilities
                for j in range(i): #j=indices of the existing nodes
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
                    hypDistances[j] = hypDist(zeta,Coords[i,:],Coords[j,:],r[i],r[j])
                    try:
                        expon = math.exp(zeta*(hypDistances[j]-R_i)/(2*T))
                    except OverflowError:
                        expon = float('inf')
                    p_i[j] = 1/(1+expon)
                targets = list(range(i)) #list of the possible target node indices
                #choose m nodes according to the probabilities in p_i[j]
                for k in range(m):
                    norm = np.sum(p_i[targets])
                    rnd = rand.random() #rndE[0,1)
                    v = 0 #iteration variable
                    w = p_i[targets[v]]/norm
                    #w=upper boundary of the interval which belongs to the currently examined event (connection to the first target)
                    #if w[v]<=rnd<w[v+1], then the (v+1)th node is chosen among the possible targets
                    #have to find the first w upper boundary, which is larger than the uniformly distributed random number rnd
                    while w<=rnd:
                        v = v+1
                        w = w+p_i[targets[v]]/norm
                    edgeList.append((targets[v],i))
                    del targets[v] #the chosen node is not a target any more
            G.add_edges_from(edgeList) #add the new edges to the graph

        #one graph has been created
        if isWeighted:
            weighting(zeta,G,Coords,r)
        print("Graph "+str(g)+" is ready.")
        Rlist.append(Rvalues)
        GraphList.append(G)
        radialCoordList.append(r)
        CartCoordList.append(Coords)
        angularGroupIDlist.append(angularGroupIDs)
    #all graphs are ready
    return [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist]



def d3nPSO_TaboveLimit(zeta,N,m,beta,T,probOfPeaks,muValues,kappaValues,numberOfGraphs,numORan,isWeighted):
    Rlist = [] #Rlist[g]=list with the cutoff distance of the connection probability in each step of the growth of the gth graph
    GraphList = [] #list of the generated NetworkX Graphs
    radialCoordList = [] #radialCoordList[g]=NumPy array with the radial coordinates of the gth graph's nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
    CartCoordList = [] #CartCoordList[g]=N-by-3 NumPy array containing the Cartesian coordinates describing the positions of the gth graph's nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
    angularGroupIDlist = []

    #Calculate the cutoff distance of the connection probability in each step of the network growth (the cutoff distance is the same for each graph)
    Rvalues = cutoffDists(zeta,N,m,beta,T,numORan) #the cutoff distance of the connection probability for node i is Rvalues[i]

    for g in range(numberOfGraphs):
        #Initialization
        G = nx.Graph() #the initial empty network
        r = np.zeros(N) #radial coordinates of the N nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2
        r_initial = (2*T*(3-1)/zeta)*np.log(np.linspace(1,N,N)) #d=3 #initial radial coordinates of the N nodes ((2*T/zeta)*ln(timestep)) in the native representation
        Coords = np.zeros((N,3)) #Cartesian coordinates of the N nodes in the native representation of the 3-dimensional hyperbolic space with curvature -zeta^2
        angularGroupIDs = -np.ones(N,dtype=int)

        #In each timestep one new node is added to the graph.
        #The time is i+1 when node i appears!!!
        #0<T: the new node i chooses its neighbors according to the probability distribution p_i[j]

        #Add the first two nodes (indexed by 0 and 1) to the graph
        G.add_node(0)
        r[0] = r_initial[0] #set the current radial coordinate of node 0 to its initial value
        Coords[0,:] = np.zeros(3) #set the initial Cartesian coordinates of node 0
        print("Node 0 does not create connections - the value of R_0 is irrelevant.")
        G.add_node(1)
        r[1] = r_initial[1] #set the current radial coordinate of node 1 to its initial value
        [angularGroupIDs[1],Coords[1,:]] = generateCartCoords(r[1],probOfPeaks,muValues,kappaValues) #generate the initial Cartesian coordinates of node 1
        #simulate popularity fading, if needed
        if beta<1: #there is popularity fading
            r[0] = beta*r_initial[0]+(1-beta)*r_initial[1]
            [angularGroupIDs[0],Coords[0,:]] = generateCartCoords(r[0],probOfPeaks,muValues,kappaValues) #update the Cartesian coordinates of node 0
            #Note that initially the length of node 0's coordinate array is 0 (r[0]=0), therefore the below applied update rule Coords[0,:] = r_0_new*Coords[0,:]/r[0] would not work here.
        #edge creation
        G.add_edges_from([(0,1)]) #add the new edge to the graph
        print("Connect node 1 to all of the previously appeared nodes - the value of R_1 is irrelevant.")

        #Continue the network growth
        for i in range(2,N):
            print("Node "+str(i))
            G.add_node(i) #add the new node to the graph
            r[i] = r_initial[i] #set the current radial coordinate of node i to its initial value
            [angularGroupIDs[i],Coords[i,:]] = generateCartCoords(r[i],probOfPeaks,muValues,kappaValues) #generate the Cartesian coordinates of node i

            edgeList = [] #list of the edges (node identifier tuples) created at the appearance of node i
            #create m edges connecting node i to the previously appeared nodes
            #note that at time i+1 the number of previously appeared nodes is i, therefore the maximum number of new edges is i
            if i<=m: #connect the new node to all of the previously appeared nodes
                print("Connect node "+str(i)+" to all of the previously appeared nodes - the value of R_"+str(i)+" is irrelevant.")
                for j in range(i): #j=indices of the previously appeared nodes
                    edgeList.append((j,i))
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
            else: #m<i: connect the new node to m previously appeared nodes chosen according to the probability distribution p_i[j]
                R_i = Rvalues[i] #set the cutoff distance of the connection probability
                #determine the connection probabilities
                hypDistances = np.zeros(i) #for the hyperbolic distances between the new node and the already existing nodes
                p_i = np.zeros(i) #for the connection probabilities
                for j in range(i): #j=indices of the existing nodes
                    if beta<1: #there is popularity fading (note that otherwise r[0]=0, therefore the division by r[j] would be problematic for j=0)
                        r_j_new = beta*r_initial[j]+(1-beta)*r_initial[i]
                        Coords[j,:] = r_j_new*Coords[j,:]/r[j] #Coords[j,:]/r[j] is a unit vector
                        r[j] = r_j_new
                    hypDistances[j] = hypDist(zeta,Coords[i,:],Coords[j,:],r[i],r[j])
                    try:
                        expon = math.exp(zeta*(hypDistances[j]-R_i)/(2*T))
                    except OverflowError:
                        expon = float('inf')
                    p_i[j] = 1/(1+expon)
                targets = list(range(i)) #list of the possible target node indices
                #choose m nodes according to the probabilities in p_i[j]
                for k in range(m):
                    norm = np.sum(p_i[targets])
                    rnd = rand.random() #rndE[0,1)
                    v = 0 #iteration variable
                    w = p_i[targets[v]]/norm
                    #w=upper boundary of the interval which belongs to the currently examined event (connection to the first target)
                    #if w[v]<=rnd<w[v+1], then the (v+1)th node is chosen among the possible targets
                    #have to find the first w upper boundary, which is larger than the uniformly distributed random number rnd
                    while w<=rnd:
                        v = v+1
                        w = w+p_i[targets[v]]/norm
                    edgeList.append((targets[v],i))
                    del targets[v] #the chosen node is not a target any more
            G.add_edges_from(edgeList) #add the new edges to the graph

        #one graph has been created
        if isWeighted:
            weighting(zeta,G,Coords,r)
        print("Graph "+str(g)+" is ready.")
        Rlist.append(Rvalues)
        GraphList.append(G)
        radialCoordList.append(r)
        CartCoordList.append(Coords)
        angularGroupIDlist.append(angularGroupIDs)
    #all graphs are ready
    return [Rlist,GraphList,radialCoordList,CartCoordList,angularGroupIDlist]



#This function calculates the cutoff distance of the connection probability in each step of the growth of a graph with given input parameters d, zeta, N, m, beta and T.
#numORan is a string specifying what type of calculations will be used:
    #if numORan=='num', the cutoff distances will be calculated by solving the equation m=i*P(i) numerically
    #if numORan=='an', the cutoff distances will be calculated analytically, using approximating formulas
def cutoffDists(zeta,N,m,beta,T,numORan):
    d=3 #3-dimensional space
    if numORan=='num':
        matEng = matlab.engine.start_matlab()
        matEng.cd(os.getcwd())
        Rvalues_matlab = matEng.cutOffDistsFor1Graph(d,zeta,N,m,beta,T,nargout=1) #a list of lists with one element
        Rvalues = [RinList[0] for RinList in Rvalues_matlab]
        matEng.quit()
        return Rvalues #the cutoff distance of the connection probability for node i is Rvalues[i]
    elif numORan=='an':
        if d%2==0: #d is even
            eta = (math.factorial(d/2-1)*math.factorial((d-2)/2)*math.pow(2,d-2)) / (math.factorial(d-2)*math.pi)
        else:
            eta = math.factorial(d-1) / (math.factorial((d-1)/2-1)*math.factorial((d-1)/2)*math.pow(2,d-1))
        Rvalues = np.zeros(N) #the cutoff distance of the connection probability for node i is Rvalues[i]
        if beta==1/(d-1):
            if T<1/(d-1):
                for i in range(m+1,N): #The nodes with index i<=m connect to all of the previously appeared nodes, therefore the cutoff distance is irrelevant for these nodes -> set it to 0
                    #Here i=0,1,...,N-1, thus the time when node i appears is i+1!
                    Rvalues[i] = (2/(zeta*(d-1)))*math.log(m*math.sin((d-1)*T*math.pi)*math.pow(i+1,2*d-3)/(eta*math.pi*math.pow(2,d-1)*T*math.log(i+1)))
                return list(Rvalues)
            elif 1/(d-1)<T:
                if d==2:
                    numericFactor = math.pow(2,1/T)*math.pow(math.pi,1-1/T)/(1-1/T)
                else:
                    matEng = matlab.engine.start_matlab()
                    matEng.cd(os.getcwd())
                    numericFactor = matEng.numIntegralInAnalyticCutoffDist(d,T,nargout=1)
                    matEng.quit()
                for i in range(m+1,N): #The nodes with index i<=m connect to all of the previously appeared nodes, therefore the cutoff distance is irrelevant for these nodes -> set it to 0
                    #Here i=0,1,...,N-1, thus the time when node i appears is i+1!
                    Rvalues[i] = (2*T/zeta)*math.log(m*math.pow(i+1,2*d-3)/(eta*numericFactor*math.log(i+1)))
                return list(Rvalues)
            else:
                print("T cannot be negative or equal to 1/(d-1)!")
        elif 0<beta and beta<=1:
            if T<1/(d-1):
                for i in range(m+1,N): #The nodes with index i<=m connect to all of the previously appeared nodes, therefore the cutoff distance is irrelevant for these nodes -> set it to 0
                    #Here i=0,1,...,N-1, thus the time when node i appears is i+1!
                    Rvalues[i] = (2/(zeta*(d-1)))*math.log(m*math.sin((d-1)*T*math.pi)*(1-(d-1)*beta)/(eta*math.pi*math.pow(2,d-1)*T*(math.pow(i+1,3-2*d)-math.pow(i+1,(d-1)*(beta-2)))))
                return list(Rvalues)
            elif 1/(d-1)<T:
                if d==2:
                    numericFactor = math.pow(2,1/T)*math.pow(math.pi,1-1/T)/(1-1/T)
                else:
                    matEng = matlab.engine.start_matlab()
                    matEng.cd(os.getcwd())
                    numericFactor = matEng.numIntegralInAnalyticCutoffDist(d,T,nargout=1)
                    matEng.quit()
                for i in range(m+1,N): #The nodes with index i<=m connect to all of the previously appeared nodes, therefore the cutoff distance is irrelevant for these nodes -> set it to 0
                    #Here i=0,1,...,N-1, thus the time when node i appears is i+1!
                    Rvalues[i] = (2*T/zeta)*math.log(m*(1-(d-1)*beta)/(eta*numericFactor*(math.pow(i+1,3-2*d)-math.pow(i+1,(d-1)*(beta-2)))))
                return list(Rvalues)
            else:
                print("T cannot be negative or equal to 1/(d-1)!")
        else:
            print('Error: The popularity fading parameter beta must fall in the range (0,1]!')
    else:
        print('Set numORan to \'num\' or \'an\'!')



#This function generates the Cartesian coordinates of a node in the native representation of the 3-dimensional hyperbolic space with the directions sampled from numOfPeaks number of von Mises–Fisher distributions. The radial coordinate of the given node in the native representation will be set to r.
#probOfPeaks is a NumPy array with numOfPeaks number of elements, each of which equal to the probability that the angular coordinates of a node will be sampled from the distribution corresponding to a given peak; the sum of the elements must be 1
#muValues is a list consisting of numOfPeaks number of lists of 3 elements, where each list of 3 elements corresponds to the Cartesian unit vector of the mean direction of a given peak in the angular distribution
#kappaValues is a list with numOfPeaks number of real positive elements, each of which corresponding to the concentration of the given peak (the greater the value of kappa, the higher the concentration of the distribution around the mean direction mu; a smaller kappa leads to a broader angular patch; kappa=0 pertains to the uniform angular distribution of the nodes)
def generateCartCoords(r,probOfPeaks,muValues,kappaValues):
    #sample the peak
    rnd = rand.random() #rnd is sampled uniformly random from the interval [0,1)
    v = 0 #iteration variable
    w=probOfPeaks[v]
    #w=upper boundary of the interval which belongs to the currently examined event (choosing the distribution with parameters muValues[v] and kappaValues[v])
    #if w[v]<=rnd<w[v+1], then the chosen distribution's (peak's) index is v+1
    #have to find the first w upper boundary, which is larger than the uniformly distributed random number rnd
    while w<=rnd:
        v=v+1
        w=w+probOfPeaks[v]
    muSelected=muValues[v]
    kappaSelected=kappaValues[v]
    vmfInstance = dist.VonMisesFisher(mean_direction=muSelected, concentration=kappaSelected)
    #sample a position on the unit sphere from the given distribution
    coordArray = vmfInstance.sample() #generate a single sample
    coordArray = coordArray.numpy()
    coordArrayWithLength_r = r*(coordArray/np.linalg.norm(coordArray)) #np.linalg.norm(coordArray) is not always exactly 1
    return [v,coordArrayWithLength_r]



#This function calculates the hyperbolic distance h between two nodes of the network. The NumPy vectors coords1 and coords2 both contain d Cartesian coordinates, which describe the position of the given nodes in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2. r1 and r2 denote the radial coordinates of the two nodes in the native representation, i.e. the Euclidean length of the vectors coords1 and coords2.
def hypDist(zeta,coords1,coords2,r1,r2):
    if r1==0:
        h = r2
    elif r2==0:
        h = r1
    else:
        cos_angle = np.inner(coords1,coords2)/(r1*r2) #cosine of the angular distance between the two nodes
        if cos_angle==1: #the vectors coords1 and coords2 point in the same direction; in this case the hyperbolic distance between the two nodes is acosh(cosh(r1-r2))=|r1-r2|
            h = math.fabs(r1-r2)
        elif cos_angle==-1: #the vectors coords1 and coords2 point in the opposite direction
            h = r1+r2
        else:
            argument_of_acosh = math.cosh(zeta*r1)*math.cosh(zeta*r2)-math.sinh(zeta*r1)*math.sinh(zeta*r2)*cos_angle
            if argument_of_acosh<1: #a rounding error occurred, because the hyperbolic distance h is close to zero
                print("The argument of acosh is "+str(argument_of_acosh)+", less than 1.")
                h = 0 #acosh(1)=0
            else:
                h = math.acosh(argument_of_acosh)/zeta
    return h




#This function weights the edges of a graph according to the hyperbolic distances between its nodes.
#zeta=sqrt(-K), where K<0 is the curvature of the hyperbolic space
#G is the NetworkX Graph to be weighted
#Coords is an N-by-d NumPy array. The ith row contains the Cartesian coordinates of node i in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2.
#nativeRadialCoords is a NumPy array with N elements, its jth element is the radial coordinate of node j in the native representation of the d-dimensional hyperbolic space with curvature -zeta^2.
def weighting(zeta,G,Coords,nativeRadialCoords):
    for (i,j) in G.edges():
        hypDist_ij = hypDist(zeta,Coords[i,:],Coords[j,:],nativeRadialCoords[i],nativeRadialCoords[j]) #calculate the hyperbolic distance between node i and j
        G[i][j]['weight'] = 1/(1+hypDist_ij)

