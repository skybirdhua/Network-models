# This file is used to model the directed scale-free network, in the paper of "Directed scale-free graphs" by Bela Bollobas 2003.

#TEST, whether it works in windows. IT works in both systems.

import random
import numpy as np
import matplotlib.pyplot as plt
from igraph import *
from scipy.stats.stats import pearsonr
from math import exp


def Initialize():
	outdeg = [1, 1, 1, 1]  # the list to recod the out-degree of nodes, from node 1,2, 3, ...
	indeg = [1, 1, 1, 1]  # the list to record the in-degree of nodes, from node 1, node 2
	links = [(0,1),(1,2),(2,3),(3,0)] # start from 0, to be consistent with the index of outdeg and indeg
	# all the elements sum up to the total edges
	stepRecord = ['0']*len(indeg) # record which step A/B/C is processed, 
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

def DirectedSF(t, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while count < t:
		#print '************   STEP : *****     ', count+1
		#print '******** The number of nodes is:  ', vertex
		prop_in = [0] * vertex
		sump_in = [0] * vertex
		prop_out = [0] * vertex
		sump_out = [0] * vertex

		sumIn = sumOut = 0  # temporary sum of all probabilities
		v = w = -1  # global variable

		for i in range(vertex):  # calculate proportional indegree
			prop_in[i] = (indeg[i] + theta_in) / (count + vertex + theta_in * vertex)
			#print 'Indegree proportional p of ', i, ' is: ',prop_in[i]
			sumIn = sumIn + prop_in[i]
			sump_in[i] = sumIn
			prop_out[i] = (outdeg[i] + theta_out) / (count + vertex + theta_out * vertex)
			#print 'Outdegree proportional p of ', i, ' is: ',prop_out[i]
			sumOut = sumOut + prop_out[i]
			sump_out[i] = sumOut
		#print 'Indegree proportion is ',prop_in, sump_in
		#print 'Outdegree proportion is ',prop_out, sump_out
		#random.seed(10)

		toss = random.random()  # produce a random number between 0 and 1

		if toss < alpha:
			# execute (A)
			# choose the old vertex w, proportial to indegree, (how to proportional ?)
			# print 'Step (A): ', toss
			#random.seed(1)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_in[i]):  # find the first one that are bigger than the probability
					break
			#print 'The target node is: ', i+1
			w = i
			indeg[w] += 1
			# add one node at the end of the lists
			indeg.append(0)
			outdeg.append(1)
			links.append((vertex,w))
			stepRecord.append('A')

			#print 'The edge is: ', vertex+1,'-->',w+1
			vertex += 1

		elif toss < alpha + beta:
			# add edge between old vertex v and w
			# print 'Step (B): ', toss
			#random.seed(2)
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			#print 'The source node is: ', i+1
			v = i

			#random.seed(3)
			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_in[i]:
					break
			#print 'The target node is: ', i+1
			w = i

			outdeg[v] += 1
			indeg[w] += 1
			links.append((v,w))
			stepRecord.append('B')
			#print 'The edge is: ', v+1,'-->',w+1

		else:
			# add a new vertex v and an input link
			#print 'Step (C): ', toss
			#random.seed(4)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_out[i]):  # find the first one that are bigger than the probability
					break
			#print 'The source node is: ', i+1
			v = i
			# w = random.randint(1,vertex)
			outdeg[v] += 1
			# add one node at the end of the lists
			indeg.append(1)
			outdeg.append(0)
			links.append((v,vertex))
			stepRecord.append('C')
			#print 'The edge is: ', v+1,'-->',vertex+1
			vertex += 1

		count += 1  # one more time


	#print 'In Initial Model: the outdegree sequence is, ', outdeg
	#print 'In Initial Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Initial Model: the fianl links are, ', links
	return links

	#indeg_hist,indeg_bin = np.histogram(indeg)
	#outdeg_hist,outdeg_bin = np.histogram(outdeg)
	#print 'The histogram of indegree is, ',indeg_hist, indeg_bin
	#print 'The histogram of outdegree is, ',outdeg_hist, outdeg_bin
	# the first value is correlation coefficient, the second one is p-values
	#print 'The correlation of in-/out-degree is, ',pearsonr(indeg,outdeg)

def DirectedSF_Neighbor_1(tt, vertex, alpha, beta, gama, theta_in, theta_out):
	# Main change 1: choose the source v according to the proportional of its out-degree;
	# Main change 2: choose the w from the neighbors of v (if any), proportional to the out-degree of neighbors
	#  if v has no neighbors, choose w globally, proportional to its in-degree

	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while (count < tt): #Can not be writen as t? why???
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
			#********** Main change 1: choose the source v according to the proportional of its out-degree**************
			#print '**********STEP B, ********'
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			v = i
			outdeg[v] += 1
			#print 'The source node is, ', v
			#********** Main change 2: choose the w from the neighbors of v *********
			#In the neighbors of v, choose the target w, proportional to its out-degree 
			
			#******Here, to improve and speed up!!!!************
			neighbors = list(set([t[0] for t in links if t[1]==v] + [t[1] for t in links if t[0]==v])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]
			
			prob2 = random.random() 

			# what if this node has only one neighbor? 
			if len(neighbors):
				#print 'The neighbors of ', v, ' include ,', neighbors
				outDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					outDegNeighbors[index] = outdeg[neighbor]
				#	print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_OutNeighbor = ToProbability(outDegNeighbors+[theta_out]*len(neighbors)) #what if the outdeg of the neighbor is zero???
				#print 'The neighbors ToProbability are, ', sump_OutNeighbor

				for i in range(len(neighbors)):
					if prob2 <= sump_OutNeighbor[i]:
						break
				w = neighbors[i]
			else:
				for i in range(vertex):
					if prob2 <= sump_in[i]:
						break
				w = i	
			indeg[w] += 1
			#print 'The target node is, ', w
			links.append((v,w))
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

	#print 'In  Neighbor_1 Model: the outdegree sequence is, ', outdeg
	#print 'In  Neighbor_1 Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Neighbor_1 Model: the fianl links are, ', links
	return links

def DirectedSF_Neighbor_1_phi(tt, vertex, alpha, beta, gama, theta_in, theta_out,phi):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while (count < tt): #Can not be writen as t? why???
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
			#**********Main change 1: choose the source v according to the proportional of its out-degree**************
			#print '**********STEP B, ********'
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			v = i
			outdeg[v] += 1
			#print 'The source node is, ', v

			#****** Main Change 2, a random variable to decide whether to choose the target globally or from the neighbors
			# if RV < phi, locally; otherwise globally
			#**********Main change 3: if locally, choose the w from the neighbors of node v, proportional to its out-degree 

			prob_phi = random.random() 
			#******Here, to improve and speed up!!!!************
			neighbors = list(set([t[0] for t in links if t[1]==v] + [t[1] for t in links if t[0]==v])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]
			# what if this node has only one neighbor? 
			
			if len(neighbors) and prob_phi <= phi:
				#print 'The neighbors of ', v, ' include ,', neighbors
				outDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					outDegNeighbors[index] = outdeg[neighbor]
				#	print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_OutNeighbor = ToProbability(outDegNeighbors+[theta_out]*len(neighbors)) #what if the outdeg of the neighbor is zero???
				#print 'The neighbors ToProbability are, ', sump_OutNeighbor

				prob2 = random.random()
				for i in range(len(neighbors)):
					if prob2 <= sump_OutNeighbor[i]:
						break
				w = neighbors[i]
			else:
				prob2 = random.random()
				for i in range(vertex):
					if prob2 <= sump_in[i]:
						break
				w = i	
			indeg[w] += 1
			#print 'The target node is, ', w
			links.append((v,w))
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

	#print 'In Neighbor_1_phi Model: the outdegree sequence is, ', outdeg
	#print 'In Neighbor_1_phi Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Neighbor_1_phi Model: the fianl links are, ', links
	return links

def DirectedSF_Neighbor_2(tt, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while (count < tt): #Can not be writen as t? why???
		v = w = -1  # global variable
		sump_in = ToProbability(list(map(lambda x: x[0]+x[1], zip(indeg,[theta_in]*len(indeg))))) # test OK
		sump_out = ToProbability(list(map(lambda x: x[0]+x[1], zip(outdeg,[theta_out]*len(outdeg))))) # test OK

		toss = random.random() 

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
			#**********Main change 1: choose the source v according to the proportional of its indegree**************
			#print '**********STEP B, ********'
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			v = i
			outdeg[v] += 1
			#print 'The source node is, ', v
			neighbors = list(set([t[0] for t in links if t[1]==v] + [t[1] for t in links if t[0]==v])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]

			prob2 = random.random() 

			#**********Main change 2: choose the w from the neighbors of node v *********
			#In the neighbors of node v, choose the target w, proportional to its in-degree 

			#******Here, to improve and speed up!!!!************
			# what if this node has only one neighbor? 
			if len(neighbors):
				#print 'The neighbors of ', w, ' include ,', neighbors
				inDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					inDegNeighbors[index] = indeg[neighbor] #!!!!! The main difference from Method Neighbor_1
				#	print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_InNeighbor = ToProbability(inDegNeighbors+[theta_in]*len(neighbors)) #what if the outdeg of the neighbor is zero???
				#print 'The neighbors ToProbability are, ', sump_OutNeighbor

				for i in range(len(neighbors)):
					if prob2 <= sump_InNeighbor[i]:
						break
				w = neighbors[i]
			else:
				for i in range(vertex):
					if prob2 <= sump_in[i]:
						break
				w = i	
			indeg[w] += 1
			#print 'The target node is, ', w
			links.append((v,w))
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

	#print 'In  Neighbor_2 Model: the outdegree sequence is, ', outdeg
	#print 'In  Neighbor_2 Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In  Neighbor_2 Model: the fianl links are, ', links
	return links

def DirectedSF_Neighbor_2_phi(tt, vertex, alpha, beta, gama, theta_in, theta_out,phi):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while (count < tt): #Can not be writen as t? why???
		v = w = -1  # global variable
		sump_in = ToProbability(list(map(lambda x: x[0]+x[1], zip(indeg,[theta_in]*len(indeg))))) # test OK
		sump_out = ToProbability(list(map(lambda x: x[0]+x[1], zip(outdeg,[theta_out]*len(outdeg))))) # test OK

		toss = random.random() 

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
			#**********Main change 1: choose the source v according to the proportional of its indegree**************
			#print '**********STEP B, ********'
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			v = i
			outdeg[v] += 1
			#print 'The source node is, ', v
			neighbors = list(set([t[0] for t in links if t[1]==v] + [t[1] for t in links if t[0]==v])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]
			#**********Main change 2: choose the w from the neighbors of node v *********
			#In the neighbors of node v, choose the target w, proportional to its in-degree 

			#******Here, to improve and speed up!!!!************
			# what if this node has only one neighbor? 
			prob_phi = random.random() 

			if len(neighbors) and prob_phi <= phi:
				#print 'The neighbors of ', w, ' include ,', neighbors
				inDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					inDegNeighbors[index] = indeg[neighbor] #!!!!! The main difference from Method Neighbor_1
				#	print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_InNeighbor = ToProbability(inDegNeighbors+[theta_in]*len(neighbors)) #what if the outdeg of the neighbor is zero???
				#print 'The neighbors ToProbability are, ', sump_OutNeighbor
				prob2 = random.random()

				for i in range(len(neighbors)):
					if prob2 <= sump_InNeighbor[i]:
						break
				w = neighbors[i]
			else:
				prob2 = random.random()
				for i in range(vertex):
					if prob2 <= sump_in[i]:
						break
				w = i	
			indeg[w] += 1
			#print 'The target node is, ', w
			links.append((v,w))
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

	#print 'In Neighbor_2_phi Model: the outdegree sequence is, ', outdeg
	#print 'In Neighbor_2_phi Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Neighbor_2_phi Model: the fianl links are, ', links
	return links
#Reciprocity is too low, hard to control the weight
def DirectedSF_Neighbor_Linear(ttt, vertex, alpha, beta, gama, theta_in, theta_out,lamb):
	# the lamb is used to balance between the degree and the directed links from neighbors
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while count < ttt:
		#print '************   STEP : *****     ', count+1
		#print '******** The number of nodes is:  ', vertex
		prop_in = [0] * vertex
		sump_in = [0] * vertex
		prop_out = [0] * vertex
		sump_out = [0] * vertex

		sumIn = sumOut = 0  # temporary sum of all probabilities
		v = w = -1  # global variable
		in_links = [0]*vertex # this list records the number of in-links to the selected source node
		out_links = [0]*vertex # this list records the number of out-links from the selected source node **think how to use this, to make symmetric


		for i in range(vertex):  # calculate proportional indegree
			prop_in[i] = (indeg[i] + theta_in) / (count + vertex + theta_in * vertex)
			#print 'Indegree proportional p of ', i, ' is: ',prop_in[i]
			sumIn = sumIn + prop_in[i]
			sump_in[i] = sumIn
			prop_out[i] = (outdeg[i] + theta_out) / (count + vertex + theta_out * vertex)
			#print 'Outdegree proportional p of ', i, ' is: ',prop_out[i]
			sumOut = sumOut + prop_out[i]
			sump_out[i] = sumOut
		#print 'Indegree proportion is ',prop_in, sump_in
		#print 'Outdegree proportion is ',prop_out, sump_out
		#random.seed(10)

		toss = random.random()  # produce a random number between 0 and 1

		if toss < alpha:
			# execute (A)
			# choose the old vertex w, proportial to indegree, (how to proportional ?)
			# print 'Step (A): ', toss
			#random.seed(1)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_in[i]):  # find the first one that are bigger than the probability
					break
			#print 'The target node is: ', i+1
			w = i
			indeg[w] += 1
			# add one node at the end of the lists
			indeg.append(0)
			outdeg.append(1)
			links.append((vertex,w))
			stepRecord.append('A')

			#print 'The edge is: ', vertex+1,'-->',w+1
			vertex += 1

		elif toss < alpha + beta:
			# add edge between old vertex v and w
			# print 'Step (B): ', toss
			#random.seed(2)
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			#print 'The source node is: ', i
			v = i

			#***************The main change is, having decided on the source node, to choose the target node: 
			#**************** impove the probability of the neighbors of the target node if any
			#random.seed(3)

			neighbors_list = list([t[0] for t in links if t[1]==v])
			neighbors = set(neighbors_list) # remove the multiple ones
			

			for neighbor in neighbors:
				in_links[neighbor] = neighbors_list.count(neighbor)
			
			#print 'The neighbors of %s : %s ' %(i, neighbors)
			#print 'The in-degree sequence is: ', indeg
			#print 'The in-links sequence is: ', in_links
			indeg0 = list(map(lambda x: (x[0]+theta_in)*(x[1]+lamb), zip(indeg,in_links)))
			#indeg0 = list(map(lambda x: lamb*x[0]+(1-lamb)*x[1], zip(indeg,in_links)))
			#********This is the key difference!!!!!!
			#sump_inNeighbor = ToProbability(list(map(lambda x: lamb*x[0]+(1-lamb)*x[1], zip(indeg,in_links))))
			#print 'The manipulated degree is, ', indeg0 
			sump_inNeighbor = ToProbability(indeg0)

			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_inNeighbor[i]:
					break
			#print 'The target node is: ', i
			w = i

			outdeg[v] += 1
			indeg[w] += 1
			links.append((v,w))
			stepRecord.append('B')
			#print 'The edge is: ', v+1,'-->',w+1

		else:
			# add a new vertex v and an input link
			#print 'Step (C): ', toss
			#random.seed(4)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_out[i]):  # find the first one that are bigger than the probability
					break
			#print 'The source node is: ', i+1
			v = i
			# w = random.randint(1,vertex)
			outdeg[v] += 1
			# add one node at the end of the lists
			indeg.append(1)
			outdeg.append(0)
			links.append((v,vertex))
			stepRecord.append('C')
			#print 'The edge is: ', v+1,'-->',vertex+1
			vertex += 1

		count += 1  # one more time


	#print 'In Initial Model: the outdegree sequence is, ', outdeg
	#print 'In Initial Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Initial Model: the fianl links are, ', links
	return links

	#indeg_hist,indeg_bin = np.histogram(indeg)
	#outdeg_hist,outdeg_bin = np.histogram(outdeg)
	#print 'The histogram of indegree is, ',indeg_hist, indeg_bin
	#print 'The histogram of outdegree is, ',outdeg_hist, outdeg_bin
	# the first value is correlation coefficient, the second one is p-values
	#print 'The correlation of in-/out-degree is, ',pearsonr(indeg,outdeg)

def DirectedSF_Neighbor_Expo(ttt, vertex, alpha, beta, gama, theta_in, theta_out,lamb):
	# the lamb is used to balance between the degree and the directed links from neighbors
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while count < ttt:
		#print '************   STEP : *****     ', count+1
		#print '******** The number of nodes is:  ', vertex
		prop_in = [0] * vertex
		sump_in = [0] * vertex
		prop_out = [0] * vertex
		sump_out = [0] * vertex

		sumIn = sumOut = 0  # temporary sum of all probabilities
		v = w = -1  # global variable
		in_links = [0]*vertex # this list records the number of in-links to the selected source node
		out_links = [0]*vertex # this list records the number of out-links from the selected source node **think how to use this, to make symmetric


		for i in range(vertex):  # calculate proportional indegree
			prop_in[i] = (indeg[i] + theta_in) / (count + vertex + theta_in * vertex)
			#print 'Indegree proportional p of ', i, ' is: ',prop_in[i]
			sumIn = sumIn + prop_in[i]
			sump_in[i] = sumIn
			prop_out[i] = (outdeg[i] + theta_out) / (count + vertex + theta_out * vertex)
			#print 'Outdegree proportional p of ', i, ' is: ',prop_out[i]
			sumOut = sumOut + prop_out[i]
			sump_out[i] = sumOut
		#print 'Indegree proportion is ',prop_in, sump_in
		#print 'Outdegree proportion is ',prop_out, sump_out
		#random.seed(10)

		toss = random.random()  # produce a random number between 0 and 1

		if toss < alpha:
			# execute (A)
			# choose the old vertex w, proportial to indegree, (how to proportional ?)
			# print 'Step (A): ', toss
			#random.seed(1)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_in[i]):  # find the first one that are bigger than the probability
					break
			#print 'The target node is: ', i+1
			w = i
			indeg[w] += 1
			# add one node at the end of the lists
			indeg.append(0)
			outdeg.append(1)
			links.append((vertex,w))
			stepRecord.append('A')

			#print 'The edge is: ', vertex+1,'-->',w+1
			vertex += 1

		elif toss < alpha + beta:
			# add edge between old vertex v and w
			# print 'Step (B): ', toss
			#random.seed(2)
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]:
					break
			#print 'The source node is: ', i
			v = i

			#***************The main change is, having decided on the source node, to choose the target node: 
			#**************** impove the probability of the neighbors of the target node if any
			#random.seed(3)

			neighbors_list = list([t[0] for t in links if t[1]==v]+[t[1] for t in links if t[0]==v])
			neighbors = set(neighbors_list) # remove the multiple ones
			

			for neighbor in neighbors:
				in_links[neighbor] = neighbors_list.count(neighbor)
			
			indeg0 = list(map(lambda x: (x[0]+theta_in)*exp(x[1]*lamb), zip(indeg,in_links)))
			#print 'The neighbors of %s : %s ' %(i, neighbors)
			#print 'The in-degree sequence is: ', indeg
			#print 'The in-links sequence is: ', in_links
			#********This is the key difference!!!!!!
			sump_inNeighbor = ToProbability(indeg0)

			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_inNeighbor[i]:
					break
			#print 'The target node is: ', i
			w = i

			outdeg[v] += 1
			indeg[w] += 1
			links.append((v,w))
			stepRecord.append('B')
			#print 'The edge is: ', v+1,'-->',w+1

		else:
			# add a new vertex v and an input link
			#print 'Step (C): ', toss
			#random.seed(4)
			prob = random.random()
			for i in range(vertex):
				if (prob <= sump_out[i]):  # find the first one that are bigger than the probability
					break
			#print 'The source node is: ', i+1
			v = i
			# w = random.randint(1,vertex)
			outdeg[v] += 1
			# add one node at the end of the lists
			indeg.append(1)
			outdeg.append(0)
			links.append((v,vertex))
			stepRecord.append('C')
			#print 'The edge is: ', v+1,'-->',vertex+1
			vertex += 1

		count += 1  # one more time


	#print 'In Initial Model: the outdegree sequence is, ', outdeg
	#print 'In Initial Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Initial Model: the fianl links are, ', links
	return links

	#indeg_hist,indeg_bin = np.histogram(indeg)
	#outdeg_hist,outdeg_bin = np.histogram(outdeg)
	#print 'The histogram of indegree is, ',indeg_hist, indeg_bin
	#print 'The histogram of outdegree is, ',outdeg_hist, outdeg_bin
	# the first value is correlation coefficient, the second one is p-values
	#print 'The correlation of in-/out-degree is, ',pearsonr(indeg,outdeg)

def DirectedSF_Inverse(t, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0

	while count < t:
		
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
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_out[i]: # v is the target
					break
			v = i
			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_in[i]: # w is the source
					break
			w = i
			#**********Main change: Inverse the proportional of in and out**************
			indeg[v] += 1 #outdeg[v] += 1
			outdeg[w] += 1 #indeg[w] += 1
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
	#print 'In Inverse Model: the outdegree sequence is, ', outdeg
	#print 'In Inverse Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Inverse Model: the fianl links are, ', links
	return links

def DirectedSF_Inverse_Neighbor(tt, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0
	while (count < tt): #Can not be writen as t? why???
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
			#print '**********STEP B, ********'
			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_in[i]:
					break
			w = i
			outdeg[w] += 1 #indeg[w] += 1
			#print 'The source node is, ', w
			#**********Main change 2: choose the v from the neighbors of w*********
			#In the neighbors of w, choose the target v, proportional to its out-degree 
			prob1 = random.random() 

			#******Here, to improve and speed up!!!!************
			neighbors = list(set([t[0] for t in links if t[1]==w] + [t[1] for t in links if t[0]==w])) # remove duplicate neighbors, as multiple edges exist #+ [t[1] for t in links if t[0]==w]
			# what if this node has only one neighbor? 
			if len(neighbors):
				#print 'The neighbors of ', w, ' include ,', neighbors
				outDegNeighbors = [0] * len(neighbors)
				for index,neighbor in enumerate(neighbors):
					outDegNeighbors[index] = outdeg[neighbor]
				#	print 'The out-degree of neighbor ', neighbor, ' is, ', outdeg[neighbor]
				sump_OutNeighbor = ToProbability(outDegNeighbors+[theta_out]*len(neighbors)) #what if the outdeg of the neighbor is zero???
				#print 'The neighbors ToProbability are, ', sump_OutNeighbor

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
			#print 'The target node is, ', v
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

	#print 'In Inverse Neighbor Model: the outdegree sequence is, ', outdeg
	#print 'In Inverse Neighbor Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Inverse Neighbor Model: the fianl links are, ', links
	return links

def DirectedSF_Inverse_Balance(t, vertex, alpha, beta, gama, theta_in, theta_out):
	indeg,outdeg,links,stepRecord = Initialize()
	count = 0

	while count < t:
		
		v = w = -1  # global variable
		sump_in = ToProbability(list(map(lambda x: x[0]+x[1], zip(indeg,[theta_in]*vertex)))) # test OK
		sump_out = ToProbability(list(map(lambda x: x[0]+x[1], zip(outdeg,[theta_out]*vertex)))) # test OK

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
			
			diffIn = [0]*len(indeg)
			diffOut = [0]*len(outdeg)

			for i in range(vertex):
				diff = indeg[i] - outdeg[i]
				if diff>0:
					diffIn[i] = diff
				else:
					diffOut[i] = -diff

			sump_diffIn = ToProbability(list(map(lambda x: x[0]+x[1], zip(diffIn,[theta_in]*vertex)))) # the balanced ones also have some probability to be linked
			sump_diffOut = ToProbability(list(map(lambda x: x[0]+x[1], zip(diffOut,[theta_out]*vertex)))) 


			#***********Main Change 1: choose target, proportional to the difference between outdegree and indegree***********
			prob1 = random.random()
			for i in range(vertex):
				if prob1 <= sump_diffOut[i]: # v is the target
					break
			v = i
			#***********Main Change 2: choose source, proportional to the difference between indegree and outdegree***********
			prob2 = random.random()
			for i in range(vertex):
				if prob2 <= sump_diffIn[i]: # w is the source
					break
			w = i
			indeg[v] += 1 #outdeg[v] += 1
			outdeg[w] += 1 #indeg[w] += 1
			links.append((w,v))
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
	#print 'In Inverse Balance Model: the outdegree sequence is, ', outdeg
	#print 'In Inverse Balance Model: the indegree sequence is, ', indeg
	#print 'The executed stepss are, ', stepRecord
	#print 'In Inverse Balance Model: the fianl links are, ', links
	return links

def GraphAnalysis(links):
	g = Graph(links,directed=True)
	v = g.vcount()
	e = g.ecount()
	numLoops = g.is_loop().count(True)
	numMulti = g.is_multiple().count(True)
	cc = g.transitivity_undirected()
	reci = g.reciprocity(mode='ratio')

	indeg = g.indegree()
	outdeg = g.outdegree()
	corr = pearsonr(indeg,outdeg)
	maxIndeg = max(indeg)
	maxOutdeg = max(outdeg)
	assor = g.assortativity_degree()

	inFit = power_law_fit(indeg)
	inFitAlpha = inFit.alpha # inFit # only run in Linux, due to the lower verion
	inFitPass = inFit.p > 0.05 # True # 

	outFit = power_law_fit(outdeg)
	outFigAlpha = outFit.alpha # outFit #  only run in linux, due to the lower version
	outFitPass = outFit.p > 0.05 # True # 

	# Add the function of bowtie structure, !!!!!!!

	#print 'Edge list is, ', g.get_edgelist()
	#print 'The correlation between in/out-degree is, ', corr
	#print 'The number of loops is, ', numLoops
	#print 'The number of multiple edges is, ', numMulti
	#print 'Clustering coefficient is, ', cc
	#print 'Graph reciprocity is, ', reci
	attriList = [v,e,cc,corr[0],reci,numLoops,numMulti,maxIndeg,maxOutdeg,assor,inFitAlpha,inFitPass,outFigAlpha,outFitPass]

	return attriList



