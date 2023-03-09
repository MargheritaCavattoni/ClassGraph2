from igraph import *
#from numba import jit, prange
import operator
import numpy as np
from igraph import *


def send_info(data, reads_info, i, to_label_dict):
    for k in range(data[i].shape[0]):
        # if the neighbour to be labelled don't have already a label
        to_label = int(data[i][k][0])
        if int(reads_info[to_label]) == 0:
            tl_weight = data[i][k][1]
            label = int(reads_info[i])
            if to_label in to_label_dict:
                to_label_dict[to_label][label] = tl_weight + to_label_dict[to_label].get(label, 0)
            else:
                to_label_dict[to_label] = {label: tl_weight}
    data[i] = np.array([])


# First version of the propagation algorithm
def lp1(max_iteration, data, reads_info):
    for v in range(max_iteration):
        to_label_dict = {}
        if v == 0:
            for i in range(len(data)):
                if int(reads_info[i]) != 0 and data[i].size > 0:
                    for k in range(data[i].shape[0]):
                        # if the neighbour to be labelled don't have already a label
                        send_info(data, reads_info, i, to_label_dict)
        else:
            for el in last_layer:
                if data[el].size > 0:
                    for k in range(data[el].shape[0]):
                        # if the neighbour to be labelled don't have already a label
                        send_info(data, reads_info, el, to_label_dict)

        if not to_label_dict:
            break
        
        last_layer = []
        for node_to_label in to_label_dict:
            v = list(to_label_dict[node_to_label].values())
            k = list(to_label_dict[node_to_label].keys())
            reads_info[node_to_label] = k[v.index(max(v))]
            last_layer.append(node_to_label)
        