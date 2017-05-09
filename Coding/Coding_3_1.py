import numpy as np
from _heapq import heappush, heappop, heapify
import matplotlib.pyplot as plt

def main():
    # build the graph
    # First two columns include nodes connected with edges
    # If graph is weighted, third column must include edge weights
    G, E, index_to_protein = buildGraph('sample_processed.csv', weighted=False)

    # Name of the PPI network. To be used in plot titles
    network_name = 'Sample'; 
        
    # all protein names
    all_proteins = []
    for i in range(len(G)):
        all_proteins.append(index_to_protein[i])
        
    # retrieve list of proteins with SNPs in Carl's genome
    f = open('Carl_Coding_SNP_Map.csv', 'r')
    lines = f.readlines()
    carl_snp_proteins = [p.strip() for p in lines]
    
    # calculate degree centrality
    degree_centrality = calculateDegreeCentrality(G, E)
    
    figure_counter = 1
    centrality_measure = 'Degree Centrality'
    # plot degree centrality (SNP and non-SNP) distributions
    plotDistributions(degree_centrality, all_proteins, carl_snp_proteins, network_name, centrality_measure, figure_counter)
    figure_counter += 2
    
    # calculate betweenness centrality
    betweenness_centrality = calculateBetweennessCentrality(G, E)

    centrality_measure = 'Betweenness Centrality'
    # plot betweenness centrality (SNP and non-SNP) distributions
    plotDistributions(betweenness_centrality, all_proteins, carl_snp_proteins, network_name, centrality_measure, figure_counter)

    plt.show()
    
# calculates betweenness centrality of each node
def calculateBetweennessCentrality(G, E):
    # a list to maintain betweenness values for each node in the graph
    betweenness_centrality_values = [0 for i in range(len(G))]
    total_number_of_shortest_paths = 0
    for source in range(len(G)):
        D, parent, number_of_shortest_paths = dijkstra(G, E, source)
        total_number_of_shortest_paths += number_of_shortest_paths
        updateBetweenness(betweenness_centrality_values, parent, source)
          
        print 'Node ' +str(source)+ ' processed. '

    # normalize betweenness centrality values as required
    normalized_betweenness_centrality_values = [float(betweenness_centrality_values[i]) / total_number_of_shortest_paths for i in range(len(betweenness_centrality_values))]

    return normalized_betweenness_centrality_values

# calculates degree centrality of nodes
def calculateDegreeCentrality(G, E):
    # degree of a node is the number of neighbors it has, ie, the length of its corresponding list in matrix G
    degree_centrality_values = [len(G[i]) for i in range(len(G))]
    
    return degree_centrality_values

# s is source node, 
# G is a graph matrix with G[i] = list of neighbors of i
# E is the edge weight matrix with E[i][j] = distance between i and its jth neighbor
def dijkstra(G, E, s):
    # initialize processed; processed[i] = 1 means node i completely processed after being visited; 0 unprocessed
    processed = [0 for i in range(len(G))]
    
    # initialize parent; parent[i] = j means i was reached through j, ie, j is parent of i in the shortest path
    parent = [-1 for i in range(len(G))]
    
    # initialize total distance; total distance to other vertices
    D = [np.Infinity for i in range(len(G))]

    # distance of source to itself is 0
    D[s] = 0

    # H is the heap of visited (but not completely processed) indices listed in order of distance
    # initialized with tuple (distance=0, source node)
    H = []
    
    # initially, all nodes are in the heap
    for node_id in range(len(D)):
        heappush(H, (D[node_id], node_id))

    while H:
        # closest_node is the closest, unvisited node in heap H and dist_to_closest is distance from source to this node
        (dist_to_closest, closest_node) = heappop(H)
        for neighbor in range(len(G[closest_node])):
            neighbor_id = G[closest_node][neighbor]
            new_distance = dist_to_closest + E[closest_node][neighbor]
            if new_distance < D[neighbor_id] and processed[neighbor_id] == 0:
                H.remove((D[neighbor_id], neighbor_id))
                heapify(H)
                D[neighbor_id] = new_distance
                parent[neighbor_id] = closest_node
                heappush(H, (new_distance, neighbor_id))

        processed[closest_node] = 1
       
    # count the number of detected shortest paths
    number_of_shortest_paths = len(parent) - parent.count(-1) + 1
     
    # return distances and parent lists
    return D, parent, number_of_shortest_paths

# reads a CSV file with edges and confidence
# returns G, E, and node_names
# G is a graph matrix (list of lists) with G[i] = list of neighbors of i
# E is the edge weight matrix (lists of lists) with E[i][j] = distance between i and its jth neighbor
# index_to_protein is a dictionary that maps index of each node to its protein name; index_to_protei[10] = "BRCA1" means node 10 corresponds to BRCA1 protein
def buildGraph(filename, weighted):
    G = []; E = []
    index_to_protein = {}
    
    # used to speed up building the graph
    protein_to_index = {}
    
    protein_counter = 0
    edge_counter = 0 
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            elements = line.strip().split(',')
            u = elements[0]; v = elements[1]; 
            edge_counter += 1
            
            if weighted:
                # For weighted networks, each edges weight should be in the 3rd column
                weight = elements[2]; 
            else:
                # For unweighted networks, edge weight = 1; 
                weight = 1

            # check if the protein have been encountered before. If not, create new nodes for them in the graph.
            index_u = protein_to_index.get(u)
            if index_u is None:
                index_u = protein_counter
                protein_to_index[u] = index_u
                index_to_protein[index_u] = u
                G.append([]); E.append([])
                protein_counter += 1

            index_v = protein_to_index.get(v)
            if index_v is None:
                index_v = protein_counter
                protein_to_index[v] = index_v
                index_to_protein[index_v] = v
                G.append([]); E.append([])
                protein_counter += 1
            
            # add edges to the graph
            G[index_u].append(index_v); G[index_v].append(index_u)
            E[index_u].append(weight); E[index_v].append(weight)

    print 'Graph constructed | ' +str(len(G))+ ' nodes, ' +str(edge_counter)+ ' edges.'
    return G, E, index_to_protein

# updates betweenness values after running Dijkstra and tracking the resulting parent list
def updateBetweenness(betweenness, parent, source):
    for current_node in range(source+1, len(betweenness)):
        # build the path back to source; stops when it reaches source and parent value = -1
        p = parent[current_node]
        
        if p != -1:
            betweenness[current_node] += 1
            
        while p != -1:
            betweenness[p] += 1
            p = parent[p]

# plot the distributions
def plotDistributions(betweenness, all_proteins, snp_proteins, network_name, graph_measure, figure_counter):
    snp_distribution = []
    nonsnp_distribution = []
    for i in range(len(all_proteins)):
        if all_proteins[i] in snp_proteins:
            snp_distribution.append(betweenness[i])
        else:
            nonsnp_distribution.append(betweenness[i])

    plt.figure(figure_counter)
    plt.xlabel(graph_measure, fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(str(graph_measure)+ ' Distribution \n SNP Proteins | ' +str(network_name)+ ' Network', fontsize=15, fontweight='bold')
    n, bins, patches = plt.hist(snp_distribution, 50, facecolor='green')
    
    plt.figure(figure_counter+1)
    plt.xlabel(graph_measure, fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(str(graph_measure)+ ' Distribution \n Non-SNP Proteins | ' +str(network_name)+ ' Network', fontsize=15, fontweight='bold')
    n, bins, patches= plt.hist(nonsnp_distribution, 50, facecolor='green')

main()