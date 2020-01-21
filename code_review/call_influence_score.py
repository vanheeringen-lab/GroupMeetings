# Python imports

import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from scipy.stats import mannwhitneyu


G = read_network("net_lung.txt")

for (u, v, d) in G.edges(data=True):
    print(u, v, d)

# EGR1 ENO1 {'weight': 0.9785623969009984, 'n': 1}
# EGR1 ERRFI1 {'weight': 0.946418743631742, 'n': 1}
# EGR1 SSU72 {'weight': 0.9418193678361, 'n': 1}
# EGR1 GNB1 {'weight': 0.9230028423551236, 'n': 1}
# EGR1 SKI {'weight': 0.908861915755856, 'n': 1}
# EGR1 AURKAIP1 {'weight': 0.8948278274298106, 'n': 1}
# EGR1 UBE2J2 {'weight': 0.8838531653496489, 'n': 1}



# expression = {'score': {'TSPAN6': 0.5953761091000407,
#                         'TNMD': 0.3838428776752495,
#                         'DPM1': 0.03578457513997305,
#                         'SCYL3': 0.30432155119404164},
#             'fc': {'TSPAN6': 0.5953761091000407,
#                         'TNMD': 0.3838428776752495,
#                         'DPM1': 0.03578457513997305,
#                         'SCYL3': 0.30432155119404164}}


def targetScore(node, G, max_degree=3, expression=None):
    """Calculate the target score."""

    if expression is None:
        expression = {"score":{}, "fc":{}}

    total_score = 0

    # get all targets up to max_degree degree 
    lengths, paths = nx.single_source_dijkstra(G, node, weight='n')
    targets = [t for t in lengths if 0 < lengths[t] <= max_degree]

    # get shortest paths based on edge weight 
    lengths, paths = nx.single_source_dijkstra(G, node, weight='weight')  
    
    # calculate target score
    for target in targets:
        path = paths[target]

        # outdegree of parent node of the target
        d = np.log(G.out_degree(path[-2]) + 1)
        # d = G.out_degree(path[-2])    

        # the level (or the number of steps) that gene is away from transcription factor
        l = len(path)   

        # expression score of the target
        g = expression["score"].get(target, 0)  
        
        # Weight is cumulative product of probabilities
        weight = [G[s][t]['weight'] for s,t in zip(path[:-1], path[1:])]

        # cumulative sum of weight 
        weight =  np.cumprod(weight)[-1]

        # score = g / l / d * weight
        score = g / l * weight
        total_score += score

    # Get Mann-Whitney U p-value of direct targets vs. non-direct targets
    direct_targets = [n for n in G[node] if n in expression["fc"]]
    non_direct_targets = [n for n in list(G.nodes) if n in expression["fc"] and n not in direct_targets]

    target_fc = [expression["fc"][t] for t in direct_targets]
    non_target_fc = [expression["fc"][t] for t in non_direct_targets]

    pval = mannwhitneyu(target_fc, non_target_fc)[1]
    target_fc_diff = np.mean(target_fc) - np.mean(non_target_fc)

    return node, total_score, G.out_degree(node), len(targets),  expression["fc"].get(node, 0), pval, target_fc_diff
    # factor, targetScore, directTargets, totalTargets, Gscore, pval, target_fc


