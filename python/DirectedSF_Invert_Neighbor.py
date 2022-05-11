import random
import numpy as np
import matplotlib.pyplot as plt
from igraph import *
from scipy.stats.stats import pearsonr



def Initialize():
	outdeg = [1, 1, 1, 1]  # the list to recod the out-degree of nodes, from node 1,2, 3, ...
	indeg = [1, 1, 1, 1]  # the list to record the in-degree of nodes, from node 1, node 2
	links = [(0,1),(1,2),(2,3),(3,1)] # start from 0, to be consistent with the index of outdeg and indeg
	# all the elements sum up to the total edges
	stepRecord = ['0']*t0 # record which step A/B/C is processed, 
	return indeg,outdeg,links,stepRecord

def ToProbability(list):
	#input a list of numbers,
	#output the probability of each numbers	
	prob = [0] * len(list)
	sumProb = [0] * len(list)
	sumList = 0

	for i,element in enumerate(list):
		sumList += element

	sumT = 0

	for i, element in enumerate(list):
		prob[i] = element/float(sumList)
		sumT += prob[i]
		sumProb[i] = sumT

	return sumProb

def DirectedSF_Inverse_Neighbor(tt, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	print 'In DirectedSF_Inverse_Neighbor Method, the t is, ', tt

	while (count < tt):
		v = w = -1  # global variable
		sump_in = ToProbability(list(map(lambda x: x[0]+x[1], zip(indeg,[theta_in]*len(indeg))))) # test OK
		sump_out = ToProbability(list(map(lambda x: x[0]+x[1], zip(outdeg,[theta_out]*len(outdeg))))) # test OK

		toss = random.random()  # produce a random number between 0 and 1

		if toss < alpha:
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_in[i]):  # find the first one that are bigger than the probability
					break
			w = i
			indeg[w] += 1
			indeg.append(0)
			outdeg.append(1)
			links.append((vertex,w))
			stepRecord.append('A')
			vertex += 1

		elif toss < alpha + beta:
			#**********Main change 1: choose the source w according to the proportional of its indegree**************
			print '**********STEP B, ********'
			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_in[i]:
					break
			w = i
			outdeg[w] += 1 #indeg[w] += 1
			print 'The source node is, ', w
			#**********Main change 2: choose the v from the neighbors of w*********
			#In the neighbors of w, choose the target v, proportional to its out-degree 
			prob1 = random.random() 
			neighbors = list(set([t[0] for t in links if t[1]==w] + [t[1] for t in links if t[0]==w])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]
			# what if this node has only one neighbor? 
			if len(neighbors):
				print 'The neighbors of ', w, ' include ,', neighbors
				outDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					outDegNeighbors[index] = outdeg[neighbor]
					print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_OutNeighbor = ToProbability(outDegNeighbors)
				print 'The neighbors ToProbability are, ', sump_OutNeighbor

				for i in range(len(neighbors)):
					if prob1 <= sump_OutNeighbor[i]:
						break
				v = neighbors[i]
			else:
				for i in range(vertex):
					if prob1 <= sump_out[i]:
						break
				v = i	
			indeg[v] += 1 #outdeg[v] += 1
			print 'The target node is, ', v
			links.append((w,v)) #links.append((v+1,w+1))
			stepRecord.append('B')

		else:
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_out[i]):  # find the first one that are bigger than the probability
					break
			v = i
			outdeg[v] += 1
			indeg.append(1)
			outdeg.append(0)
			links.append((v,vertex))
			stepRecord.append('C')
			vertex += 1
		count += 1  # one more time

	print 'The outdegree sequence is, ', outdeg
	print 'The indegree sequence is, ', indeg
	print 'The executed stepss are, ', stepRecord
	print 'The links are, ', links
	#return links

t = 10  # the total number of time steps, from small to big, 6, 10, 20, 100, 1000, 10000
t0 = 4  # the initial edges in G0

alpha = 0.3  # random.random() # the probability to execute (A)
beta = 0.4  # random.random() # the probability to execute (B)
gama = 0.3  # random.random() # the probability to execute (C)
# alpha + beta + gama = 1

theta_in = 0.1  # the parameter that controls the popularity of indegree
theta_out = 0.0

vertex = 4  # the initial number of nodes, can be changed

DirectedSF_Inverse_Neighbor(t,vertex,alpha,beta,gama,theta_in,theta_out)
