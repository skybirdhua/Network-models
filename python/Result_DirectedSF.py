# This file is used to generate the averaged results after a lot of iterations, by running the methods defined in Directed_scale_free.py

# All the methods defined in that file are simple variants based on the original model of Directed scale-free, descriped in the paper

from Directed_scale_free import *
#Initialize,ToProbability,DirectedSF,DirectedSF_Neighbor_1,DirectedSF_Neighbor_1_phi,DirectedSF_Neighbor_2,DirectedSF_Neighbor_2_phi,DirectedSF_Inverse,DirectedSF_Inverse_Neighbor,DirectedSF_Inverse_Balance,DirectedSF_Neighbor_Linear,GraphAnalysis


t = 9996 #36653  # the total number of time steps, from small to big, 6, 10, 20, 100, 1000, 10000
t0 = 4  # the initial edges in G0

alpha = 0.069  # random.random() # the probability to execute (A)
beta = 0.862  # random.random() # the probability to execute (B)
gama = 0.069  # random.random() # the probability to execute (C)
# alpha + beta + gama = 1
# alpha + gama = 0.13776359221976703 in MDD, beta = 0.8622364077802329, alpha = gama = 0.06888179610988351
# alpha + gama = 0.11556837556004916 in Dep, beta = 0.8844316244399508, alpha = gama = 0.05778418778002458

#lamb = 0.0000000000000001 # This new parameter is used to balance the in-degree and in-links (to the target node)
theta_in = 0.1 #1/(alpha+beta) * (0.3)  # the mean degree 1/(alpha+beta) but consider the number of edges #4 #2 #1 #0.1 #9.007  # the parameter that controls the popularity of indegree
theta_out = 0.1 #1/(alpha+beta) * (0.3) # the mean degree 1/(alpha+beta) #3 #2 #1 #0.1 #8.177

phi = 0.7 # this parameter serve as the control parameter, decides to local or globally select the target.
#Test the most propriate phi

vertex = 4  # the initial number of nodes, can be changed


iteration = 100
nodes = [0]*iteration
edges = [t+t0]*iteration
clucs = [0]*iteration
corrs = [0]*iteration
recis = [0]*iteration
loops = [0]*iteration
multi = [0]*iteration
maxIn = [0]*iteration
maxOut = [0]*iteration
assor = [0]*iteration
inFit = [0]*iteration
inPass = [False]*iteration
outFit = [0]*iteration
outPass = [False]*iteration

#Before run the methods, remember to change the printing xxx 
print 'The current phi is ', phi
print '************0.0 Executing Original Model**************'
for i in range(iteration):
	links = DirectedSF(t,vertex,alpha,beta,gama,theta_in,theta_out)
# 	#print '************Totally : ', t+t0, 'edges and ', nodes1, ' nodes'
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration

print '************1.1. Executing Original Neighbor_1 Model**************'
for i in range(iteration):
	links = DirectedSF_Neighbor_1(t,vertex,alpha,beta,gama,theta_in,theta_out)
# 	#print '************Totally : ', t+t0, 'edges and ', nodes1, ' nodes'
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration


print '************1.2. Executing Original Neighbor_1 phi Model**************'
for i in range(iteration):
	links = DirectedSF_Neighbor_1_phi(t,vertex,alpha,beta,gama,theta_in,theta_out,phi)
# 	#print '************Totally : ', t+t0, 'edges and ', nodes1, ' nodes'
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration



print '************2.1. Executing Original Neighbor_2 Model**************'
for i in range(iteration):
	links = DirectedSF_Neighbor_2(t,vertex,alpha,beta,gama,theta_in,theta_out)
# 	#print '************Totally : ', t+t0, 'edges and ', nodes1, ' nodes'
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration

print '************2.2. Executing Original Neighbor_2_Phi Model**************'
for i in range(iteration):
	links = DirectedSF_Neighbor_2_phi(t,vertex,alpha,beta,gama,theta_in,theta_out,phi)
# 	#print '************Totally : ', t+t0, 'edges and ', nodes1, ' nodes'
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration

print '************3. Executing Inversed Neighbors Model**************'
for i in range(iteration):
	links = DirectedSF_Inverse_Neighbor(t,vertex,alpha,beta,gama,theta_in,theta_out)
	attriList = GraphAnalysis(links)
	nodes[i] = attriList[0]
	clucs[i] = attriList[2]
	corrs[i] = attriList[3]
	recis[i] = attriList[4]
	loops[i] = attriList[5]
	multi[i] = attriList[6]
	maxIn[i] = attriList[7]
	maxOut[i] = attriList[8]
	assor[i] = attriList[9]
	inFit[i] = attriList[10]
	inPass[i] = attriList[11]
	outFit[i] = attriList[12]
	outPass[i] = attriList[13]

print 'Average node number is, ', sum(nodes)/iteration
print 'Average edge number is, ', sum(edges)/iteration
print 'Average clustering coefficient is, ', sum(clucs)/iteration
print 'Average correlation is, ', sum(corrs)/iteration
print 'Average reciprocity is, ', sum(recis)/iteration
print 'Average loops is, ', sum(loops)/iteration
print 'Average mutiedge is, ', sum(multi)/iteration
print 'Average max indegree is, ', sum(maxIn)/iteration
print 'Average max outdegree is, ', sum(maxOut)/iteration
print 'Average assortativity is, ', sum(assor)/iteration
print 'Average indegree exponent is, ', sum(inFit)/iteration
print 'Indegree fit pass rate is, ', inPass.count(True)/iteration
print 'Average outdegree exponent is, ', sum(outFit)/iteration
print 'Outdegree fit pass rate is, ', outPass.count(True)/iteration

print '##########################'

# print '************4. Executing Inversed balanced Model**************'
# for i in range(iteration):
# 	links = DirectedSF_Inverse_Balance(t,vertex,alpha,beta,gama,theta_in,theta_out)
# 	attriList = GraphAnalysis(links)
# 	nodes4[i] = attriList[0]
# 	corrs4[i] = attriList[1]
# 	recis4[i] = attriList[2]
# 	loops4[i] = attriList[3]
# 	multi4[i] = attriList[4]
# 	maxIn4[i] = attriList[5]
# 	maxOut4[i] = attriList[6]
# 	assor4[i] = attriList[7]
# 	inFit4[i] = attriList[8]
# 	inPass4[i] = attriList[9]
# 	outFit4[i] = attriList[10]
# 	outPass4[i] = attriList[11]
# print 'Average node number is, ', sum(nodes4)/iteration
# print 'Average correlation is, ', sum(corrs4)/iteration
# print 'Average reciprocity is, ', sum(recis4)/iteration
# print 'Average loops is, ', sum(loops4)/iteration
# print 'Average mutiedge is, ', sum(multi4)/iteration
# print 'Average max indegree is, ', sum(maxIn4)/iteration
# print 'Average max outdegree is, ', sum(maxOut4)/iteration
# print 'Average assortativity is, ', sum(assor4)/iteration
# print 'Average indegree exponent is, ', sum(inFit4)/iteration
# print 'Indegree fit pass rate is, ', inPass4.count(True)/iteration
# print 'Average outdegree exponent is, ', sum(outFit4)/iteration
# print 'Outdegree fit pass rate is, ', outPass4.count(True)/iteration





# analyze the network and obtain the properties, including in-deg distribution, out-deg distribution,
# correlation between in-deg and out-deg, the joint in-deg/out-deg distribution
# other network properties, including, the average degree, the max/min degree, the clustering coefficient,
# the reciprocity, the Bow Tie structure, etc.
# the motif count, 
# the # of self-loops, the # of multiple edges,
# the weighted reciprocity


# To have 1000 networks, and take the average of these properties
