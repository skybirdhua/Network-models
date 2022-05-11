# This file is used to deal with the real-world empirical data,
# The empiracle data is from the two groups of douban, one is MDD and the other is Depression

# Read the communication links

import csv
from igraph import *
from Directed_scale_free import GraphAnalysis


def ToLinks(name):
	# e.g. name=group151898-new-network
	rpath = '../Network Datasets/Networks/'+name+'.csv'
	wpath = 'Links_'+name+'.txt'

	links=[]
	source=[]
	target=[]
	topicID=[]
	words=[]
	time=[]
	users=[] # record the user's ID, and use its index as the vid

	with open(rpath,'rb') as netfile:
		rows = csv.reader(netfile,delimiter=',')
		for row in rows:
			if not (row[0] in users):
				users.append(row[0])
			if not (row[1] in users):
				users.append(row[1])
			source.append(row[0])
			target.append(row[1])
			topicID.append(row[3])
			words.append(row[4])
			time.append(row[5])

	with open(wpath,'w') as linkfile:
		for s, t in zip(source[1:], target[1:]):
			linkfile.write(str(users.index(s)))
			linkfile.write(' ')
			linkfile.write(str(users.index(t)))
			linkfile.write('\n')

	netfile.close()
	linkfile.close()



#Here the programme starts
#name = 'group151898-new-network'
name = 'fly_vs_free'
# Before that, make sure execute the following, 
# ToLinks(name)
rpath = 'Links_'+name+'.txt'
links = []

with open(rpath,'r') as linkfile:
	linkLines = linkfile.read().split('\n')
	for line in linkLines:
		link = line.split(' ')
		#if not (int(link[0]) == int(link[1])): #remove self-loops here
		links.append((int(link[0]),int(link[1])))


attriList = GraphAnalysis(links)
print 'The number of nodes is,    ',attriList[0]
print 'The number of edges is,    ',attriList[1]
print 'The clustering coeffi is,  ',attriList[2]
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


linkfile.close()


#name = 'group151898-new-network'
# The number of nodes is,     5052
# The number of edges is,     36657
# The clustering coeffi is,     0.0447391951497
# The in/out correlation is,  0.759880299208
# The reciprocity is,         0.340488554085
# The number of loops is,     0
# The number of multiedge is, 14067
# The max indegree is,        451
# The max outdegree is,       942
# The assortativity is,       0.0294356021216
# The indegree exponent is,   2.35565097084
# The indegree fit pass? :    True
# The outdegree exponent is,  2.23108271762
# The outdegree fit pass? :   True


#name = 'fly_vs_free'
# The number of nodes is,     11661
# The number of edges is,     100884
# The clustering coeffi is,     0.0425808501063
# The in/out correlation is,  0.809387192749
# The reciprocity is,         0.330572408336
# The number of loops is,     0
# The number of multiedge is, 42365
# The max indegree is,        2391
# The max outdegree is,       1553
# The assortativity is,       0.339345684059
# The indegree exponent is,   2.31990580235
# The indegree fit pass? :    True
# The outdegree exponent is,  2.46650759185
# The outdegree fit pass? :   True

