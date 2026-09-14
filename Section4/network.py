#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:51:15 2025
@author: Giorgio Ricchiuti 
@website: www.grarchive.net

Description

"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def create_network(N):
        
    ws_graph=nx.complete_graph(N)
    #The network is fully connected.
    A=nx.to_numpy_array(ws_graph) # Let's save the adjancency matrix of G
    np.fill_diagonal(A, 1) # Let's fill the diagonal
    #Let's calculate the degree for each node. In the model is crucial to calculate the effect of neighbors.  
    B=np.array(ws_graph.degree()) # array with degrees
    C=B[:,1]

    return ws_graph, A, B, C