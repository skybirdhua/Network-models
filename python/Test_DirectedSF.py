#This file is used to test new functions added in the Directed_scale_free.py 

from Directed_scale_free import *
#Initialize,ToProbability,DirectedSF,DirectedSF_Neighbor,DirectedSF_Inverse,DirectedSF_Inverse_Neighbor,DirectedSF_Inverse_Balance,DirectedSF_Neighbor_Linear,GraphAnalysis



#print Initialize() # test OK

t = 96  # the total number of time steps, from small to big, 6, 10, 20, 100, 1000, 10000
t0 = 4  # the initial edges in G0

alpha = 0.3  # random.random() # the probability to execute (A)
beta = 0.4  # random.random() # the probability to execute (B)
gama = 0.3  # random.random() # the probability to execute (C)
# alpha + beta + gama = 1

lamb = 0.1 # This new parameter is used to balance the in-degree and in-links (to the target node)
theta_in = 0.1  # the parameter that controls the popularity of indegree
theta_out = 0.1

vertex = 4  # the initial number of nodes, can be changed

print '************1. Executing Original Model**************'
links = DirectedSF(t,vertex,alpha,beta,gama,theta_in,theta_out)
attriList = GraphAnalysis(links)
print 'The number of nodes is,    ',attriList[0]
print 'The number of edges is,    ',t+t0
print 'The clustering coefficient is, ',attriList[2]
print 'The in/out correlation is, ',attriList[3]
print 'The reciprocity is,        ',attriList[4]
print 'The number of loops is,    ',attriList[5]
print 'The number of multiedge is,',attriList[6]
print 'The max indegree is,       ',attriList[7]
print 'The max outdegree is,      ',attriList[8]
print 'The assortativity is,      ',attriList[9]
print 'The indegree exponent is,  ',attriList[10]
print 'The indegree fit pass? :   ',attriList[11]
print 'The outdegree exponent is, ',attriList[12]
print 'The outdegree fit pass? :  ',attriList[13]


print '************2. Executing Original Neighbor Model**************'
links = DirectedSF_Neighbor(t,vertex,alpha,beta,gama,theta_in,theta_out)
attriList = GraphAnalysis(links)
print 'The number of nodes is,    ',attriList[0]
print 'The number of edges is,    ',t+t0
print 'The clustering coefficient is, ',attriList[2]
print 'The in/out correlation is, ',attriList[3]
print 'The reciprocity is,        ',attriList[4]
print 'The number of loops is,    ',attriList[5]
print 'The number of multiedge is,',attriList[6]
print 'The max indegree is,       ',attriList[7]
print 'The max outdegree is,      ',attriList[8]
print 'The assortativity is,      ',attriList[9]
print 'The indegree exponent is,  ',attriList[10]
print 'The indegree fit pass? :   ',attriList[11]
print 'The outdegree exponent is, ',attriList[12]
print 'The outdegree fit pass? :  ',attriList[13]


# print '************3. Executing Inversed Model**************'
# links = DirectedSF_Inverse(t,vertex,alpha,beta,gama,theta_in,theta_out)
# attriList = GraphAnalysis(links)
# print 'The number of nodes is,    ',attriList[0]
# print 'The in/out correlation is, ',attriList[1]
# print 'The reciprocity is,        ',attriList[2]
# print 'The number of loops is,    ',attriList[3]
# print 'The number of multiedge is,',attriList[4]
# print 'The max indegree is,       ',attriList[5]
# print 'The max outdegree is,      ',attriList[6]
# print 'The assortativity is,      ',attriList[7]
# print 'The indegree exponent is,  ',attriList[8]
# print 'The indegree fit pass? :   ',attriList[9]
# print 'The outdegree exponent is, ',attriList[10]
# print 'The outdegree fit pass? :  ',attriList[11]

print '************4. Executing Inversed Neighbors Model**************'
links = DirectedSF_Inverse_Neighbor(t,vertex,alpha,beta,gama,theta_in,theta_out)
attriList = GraphAnalysis(links)
print 'The number of nodes is,    ',attriList[0]
print 'The number of edges is,    ',t+t0
print 'The clustering coefficient is, ',attriList[2]
print 'The in/out correlation is, ',attriList[3]
print 'The reciprocity is,        ',attriList[4]
print 'The number of loops is,    ',attriList[5]
print 'The number of multiedge is,',attriList[6]
print 'The max indegree is,       ',attriList[7]
print 'The max outdegree is,      ',attriList[8]
print 'The assortativity is,      ',attriList[9]
print 'The indegree exponent is,  ',attriList[10]
print 'The indegree fit pass? :   ',attriList[11]
print 'The outdegree exponent is, ',attriList[12]
print 'The outdegree fit pass? :  ',attriList[13]


# print '************5. Executing Inversed balanced Model**************'
# links = DirectedSF_Inverse_Balance(t,vertex,alpha,beta,gama,theta_in,theta_out)
# attriList = GraphAnalysis(links)
# print 'The number of nodes is,    ',attriList[0]
# print 'The number of edges is,    ',t+t0
# print 'The clustering coefficient is, ',attriList[2]
# print 'The in/out correlation is, ',attriList[3]
# print 'The reciprocity is,        ',attriList[4]
# print 'The number of loops is,    ',attriList[5]
# print 'The number of multiedge is,',attriList[6]
# print 'The max indegree is,       ',attriList[7]
# print 'The max outdegree is,      ',attriList[8]
# print 'The assortativity is,      ',attriList[9]
# print 'The indegree exponent is,  ',attriList[10]
# print 'The indegree fit pass? :   ',attriList[11]
# print 'The outdegree exponent is, ',attriList[12]
# print 'The outdegree fit pass? :  ',attriList[13]

