from igraph import *
from numba import njit
import operator
import numpy as np
from igraph import *

# First version of the propagation algorithm
def lp1(max_iteration, data, reads_info):
    for v in range(max_iteration):
        # print("iter" + str(v))

        # Nodes at level i-1 send info to their neighbours at level i
        tolabel_count = 0
        for i in range(len(data)):
            if int(reads_info[i]) != 0 and len(data[i]) > 0:
                for k in range(len(data[i])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][k][0])
                    if reads_info[to_label] == 0:
                        tolabel_count += 1
                        tmp = []
                        tl_weight = data[i][k][1]
                        label = reads_info[i]
                        tmp.append(tl_weight)
                        tmp.append(label)
                        data[to_label].append(tmp)
                data[i] = []

        #print(tolabel_count)

        if tolabel_count == 0:
            break
        # For each node at level i the final label is decided

        for i in range(len(data)):
            len_line = len(data[i])
            if len_line > 1:
                possible_labels = []
                for k in range(1, len_line):
                    possible_labels.append(data[i][k])
                possible_labels = sorted(possible_labels, key=operator.itemgetter(1))
                summing_list = []
                for j in range(len(possible_labels)):
                    if len(summing_list) == 0:
                        summing_list.append(possible_labels[j])
                    elif len(summing_list) > 0 and summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                        summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + \
                                                                 possible_labels[j][0]
                    else:
                        summing_list.append(possible_labels[j])
                summing_list = sorted(summing_list, key=operator.itemgetter(0))
                reads_info[i] = summing_list[len(summing_list) - 1][1]
                for k in range(1, len_line):
                    del data[i][1]



# Second version of the propagation algorithm

def lp2(max_iteration, data, beta):
    for v in range(max_iteration):
        # print("iter" + str(v))

        # Nodes at level i-1 send info to their neighbours at level i

        tolabel_count = 0
        for i in range(len(data)):
            if int(data[i][1] != 0) and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][2][k][0])
                    if data[to_label][1] == 0:
                        tolabel_count += 1
                        tmp = []
                        tl_weight = data[i][2][k][1]
                        label = data[i][1]
                        tmp.append(tl_weight)
                        tmp.append(label)
                        data[to_label].append(tmp)
                        data[to_label][3] += 1
                        # print(data[to_label][3])
                data[i][2] = []

        #print(tolabel_count)

        if tolabel_count == 0:
            break

        # Nodes at level i exchange infromations between each other

        for i in range(len(data)):
            if int(data[i][1] == 0) and len(data[i]) > 4:
                for k in range(len(data[i][2])):
                    if data[data[i][2][k][0]][1] == 0 and len(data[data[i][2][k][0]]) > 4:
                        for j in range(4, 4 + data[i][3]):
                            new_tmp = []
                            mweight = data[i][j][0] * data[i][2][k][1] * beta
                            new_tmp.append(mweight)
                            new_tmp.append(data[i][j][1])
                            data[data[i][2][k][0]].append(new_tmp)

        # For each node at level i the final label is decided

        for i in range(len(data)):
            len_line = len(data[i])
            if len_line > 4:
                possible_labels = []
                for k in range(4, len_line):
                    possible_labels.append(data[i][k])
                possible_labels = sorted(possible_labels, key=operator.itemgetter(1))
                summing_list = []
                for j in range(len(possible_labels)):
                    if len(summing_list) == 0:
                        summing_list.append(possible_labels[j])
                    elif len(summing_list) > 0 and summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                        summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + \
                                                                 possible_labels[j][0]
                    else:
                        summing_list.append(possible_labels[j])
                summing_list = sorted(summing_list, key=operator.itemgetter(0))
                data[i][1] = summing_list[len(summing_list) - 1][1]
                data[i][3] = 0
                for k in range(4, len_line):
                    del data[i][4]

