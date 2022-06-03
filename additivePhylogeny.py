import sys
import numpy as np
import math
import networkx as nx
import matplotlib.pyplot as plt
import pylab

G = nx.Graph()
global showDataBool
global showMatrixBool

def fromInput():
    data = []
    try:
        size = float(input("Size of matrix: "))
        for i in range(size):
            row = list(map(float, input(fr"Row {i}: ").strip().split()))
            numbers = np.array(row)
            if numbers.size != size:
                quit()
            data.append(numbers)
    except:
        print("Error")
        quit()
    return np.array(data)

def fromFile(path):
    try:
        data = []
        dataFile = open(path, "r")
        row = list(map(float, dataFile.readline().strip().split()))
        size = len(row)
        data.append(np.array(row))
        i = 1
        while size > i:
            numbers = np.array(list(map(float, dataFile.readline().strip().split())))
            if numbers.size != size:
                quit()
            data.append(numbers)
            i += 1
        dataFile.close()
    except:
        print("Error")
        quit()
    return np.array(data)

def get_edge_attributes(G, name):
    edges = G.edges(data=True)
    return dict( (x[:-1], x[-1][name]) for x in edges if name in x[-1] )

def showData(data):
    if showMatrixBool:
        print(data)

def limbLength(data, j):
    limb = np.inf
    iList = [i for i in range(len(data)) if i != j]
    for k in range(len(data)):
        if k == j:
            continue
        limb = min(limb, min(data[iList,j] - data[iList,k] + data[j,k]))
    return limb//2

def findLeaves(data):
    for k in range(len(data)-1):
        arr = data[k] - data[-1]
        index = np.where(arr == data[k, -1])
        if len(index[0]) > 0:
            return (index[0][0], k)
    return None

def nearestNode(edge, weight, x, i, k):
    queue = [[i]]
    visited = set([i])
    findPath = []
    while len(queue) > 0:
        path = queue.pop()
        node = path[-1]
        visited.add(node)
        if node == k:
            findPath = path
            break
        for next_node in edge[node]:
            if next_node not in visited:
                queue.append(path+[next_node])
    dist = 0
    for k in range(len(findPath)-1):
        i, j = findPath[k], findPath[k+1]
        if dist+weight[(i, j)] > x:
            return (i, j, x-dist, dist+weight[(i, j)]-x)
        dist += weight[(i, j)]

def additivePhylogeny(data, n, innerN):
    showData(data)
    if n == 2:
        edge = {}
        edge[0] = [1]
        edge[1] = [0]
        weight = {}
        weight[(0, 1)] = data[0, 1]
        weight[(1, 0)] = data[0, 1]
        return (edge, weight, innerN)
    limb = limbLength(data, n-1)
    data[:-1,-1] -= limb
    data[-1,:-1] -= limb
    i, k = findLeaves(data)
    x = data[i, -1]
    edge, weight, innerN = additivePhylogeny(data[:-1, :-1], n-1, innerN)
    iNearest, kNearest, ix, nx = nearestNode(edge, weight, x, i, k)
    newNode = iNearest
    if ix != 0:
        newNode = innerN
        innerN += 1
        edge[iNearest].remove(kNearest) 
        edge[kNearest].remove(iNearest)
        edge[iNearest].append(newNode)
        edge[kNearest].append(newNode)
        edge[newNode] = [iNearest, kNearest]
        weight[(newNode, iNearest)] = ix
        weight[(iNearest, newNode)] = ix
        weight[(newNode, kNearest)] = nx
        weight[(kNearest, newNode)] = nx
        del weight[(iNearest, kNearest)]
        del weight[(kNearest, iNearest)]
    edge[newNode].append(n-1)
    edge[n-1] = [newNode]
    weight[(n-1, newNode)] = limb
    weight[(newNode, n-1)] = limb
    return (edge, weight, innerN)

def createTree(edge, weight, n):
    leafCount = 0
    clusterCount = 0
    for i in range(len(edge)):
        if leafCount < n:
            newLeaf = "L"+str(leafCount)
            G.add_node(newLeaf)
            leafCount +=1
            continue
        newCluster = "C"+str(clusterCount)
        G.add_node(newCluster)
        clusterCount +=1
    for i in sorted(edge):
        for j in sorted(edge[i]):
            first = "L"+str(i) if i < n else "C"+str(i-leafCount)
            second = "L"+str(j) if j < n else "C"+str(j-leafCount)
            G.add_edge(first, second, weight = weight[(i, j)])

if __name__=='__main__':
    showMatrixBool = True
    showDataBool = True
    #data = fromInput()
    data = fromFile(fr"input\test5.txt")
    try:
        n = len(data)
        edge, weight, _ = additivePhylogeny(data, n, n)
        createTree(edge, weight, n)
        pos = nx.spring_layout(G, k=2.15, iterations=20)
        edge_labels=dict([((u,v,),d['weight'])
                    for u,v,d in G.edges(data=True)])
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        if showDataBool:
            print("Nodes:",G.nodes)
            print("Edges with weights:", get_edge_attributes(G, "weight"))
        nx.draw(G, pos, node_size=500,edge_cmap=plt.cm.Reds,with_labels=True)
        pylab.show()
    except:
        print("It is not a additive tree")