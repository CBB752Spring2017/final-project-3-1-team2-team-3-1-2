---
layout: page
title: CBB752 Spring 2017
tagline: Final Project
---

3.1 - Network Analysis of the Human Protein-Protein Interactome
------------------


Table of Contents
-----------------------
1. [Introduction](#introduction)
2. [Coding](#coding)
3. [Pipeline](#pipeline)
4. [Writing](#writing)
5. [Conclusions](#conclusions)
6. [References](#references)




**Contributors**
 -Writing: Amy Zhao
 -Coding: Hussein Mohsen
 -Pipeline: Dingjue Ji

### Introduction:
#### Project Summary
Human protein-protein interaction (PPI) networks provide valuable insight into the functionality of proteins beyond what is detailed in the human proteome. 

In this project, we hope to visualize the structure of the human PPI network - which was obtained from two different databases - as well as determine the distribution of the degree and betweenness centrality measurements for the proteins within said network, using both our own proposed tool as well as the bioinformatics software Cytoscape. Furthermore, we would like to characterize any statistically significant differences between the proteins containing and not containing SNPs in Carl's genome. Finally, we also perform a hierarchical analysis of the PPI network.

For project 3.1, protein-protein interaction (PPI) data were downloaded from two different databases: The Database of Interacting Proteins (DIP) and the Molecular Interaction Database (MINT). We filtered the edges where one or two of the proteins did not have the UniProtID, which constituted less than 5% of the total proteins of interest. The former database generated 4904 distinct proteins with 7,387 interactions, while the latter yielded 4584 distinct proteins with 12,655 interactions. 


### Writing:








### Coding:

#### Documentation: 
##### Summary: 
Hussein coded a Python program that calculates the betweenness centrality and degree centrality measurements for each node in the PPI network. Figures 3-6 show the distributions for the degree and betweenness centralities outputted by Hussein's program. The observations of these distributions will be discussed in the following section. 

##### Separating SNP-Containing Proteins from Non-SNP-Containing Proteins
To help separate the proteins that contain SNPs in Carl’s genome from those that do not, the proteins’ UniProt IDs were used; that is, the proteins from the data file Carl_Coding_SNP_Map.csv – obtained from Carl’s Game of Genomes blog – were mapped to the proteins from the database using the corresponding UniProtKB ID-ENSEMBL Transcript ID pairs. The total number of SNP-containing proteins is 375 in the DIP database and 306 in the MINT database. 

##### Calculation of Degree Centrality and Betweenness Centrality
###### Centrality Measure Definitions
Next, the degree centrality and betweenness centrality measurements were calculated for each node. It is first important to define these two measures. The former is equivalent to the number of the links a node has. On the other hand, the betweenness centrality is defined as the number of times a node is located on the shortest path between two other nodes. In mathematical terms, the two measures are defined as follows: 
For a graph $G := (V, E)$ and given a node $n$: 
$$C_D(n) = deg(n)$$
$$C_B(n) = \sum_{l \ne m \ne n \in V}\frac{\sigma_{lm}(n)}{\sigma_{lm}}$$
where $\sigma_{lm}$ is the total number of shortest paths between distinct nodes $l$ and $m$. 

Given an input CSV file with three columns corresponding to protein interactor A, protein interactor B, and weights (for an example, see sample_processed.csv), Hussein's code computes the degree centrality and betweenness centrality for each node within the network. The specific documentation for each line of the code from Hussein can found below: 

Script Coding_3_1.py is rife with self-explanatory comments.

Required Packages:
==================
matplotlib and numpy


Inputs:
======= 
- A file that with the edges in the PPI network: with IDs of connected proteins in the first two columns (submitted sample_processed.csv is a selected sample from DIP network with 1000 edges and 1199 nodes/proteins)
- A file with list of protein IDs that include SNPs (submitted Carl_Coding_SNP_Map.csv is a sample file). For Carl's SNP proteins, they were identified as per per ENSEMBLE IDs in coding SNPs file and mapped to their UniProt IDs.
- Default is unweighted PPI networks. Accordingly, third column is optional: if not empty, it should include the weight of the edge so that the script accommodates weighted graphs.

Output:
=======

Degree centrality and betweenness centrality distributions plotted for proteins with and without SNPs.




#### Results: 


### Pipeline:


#### Documentation:


#### Results:









#### Conclusions:








#### References:

 References can be included here or at the end of each relevant section.
 
 
