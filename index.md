---
layout: page
title: CBB752 Spring 2017
tagline: Final Project
---

Project Title
------------------


Table of Contents
-----------------------




**Contributors**
 -Writing: Amy Zhao
 -Coding: Hussein Mohsen
 -Pipeline: Dingjue Ji

### Introduction:





### Writing:








### Coding:Calculate the degree centrality and betweenness centrality of proteins containing and not containing SNPs in Carl’s genome using a PPI file.

#### Documentation: 

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



#### Results: The distributions for both DIP and MINT networks are very similar to those by Cytoscape. They are shown in the final report submitted by Amy Zhao.



### Pipeline:


#### Documentation:


#### Results:









#### Conclusions:








#### References:

 References can be included here or at the end of each relevant section.
 
 
