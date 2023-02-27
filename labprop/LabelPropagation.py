from igraph import *
from numba import njit
import operator
import numpy as np
# First version of the propagation algorithm

def lp1(max_iteration : int, data : list):
    for v in range(max_iteration):
        # print("iter" + str(v))

        # Nodes at level i-1 send info to their neighbours at level i
        tolabel_count = 0
        for i in range(len(data)):
            if int(data[i][1]) != 0 and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][2][k][0])
                    if data[to_label][1] == 0:
                        tolabel_count += 1
                        tmp = np.array([])
                        tl_weight = data[i][2][k][1]
                        label = data[i][1]
                        np.append(tmp, tl_weight)
                        np.append(tmp, label)
                        np.append(data[to_label], tmp)
                data[i][2] = np.array([])

        #print(tolabel_count)

        if tolabel_count == 0:
            break
        # For each node at level i the final label is decided

        for i in range(len(data)):
            len_line = len(data[i])
            if len_line > 3:
                possible_labels = np.array([])
                for k in range(3, len_line):
                    possible_labels.append(data[i][k])
                #possible_labels = sorted(possible_labels, key=operator.itemgetter(1))
                possible_labels = possible_labels [possible_labels [:, 1].argsort()]
                summing_list = np.array([])
                for j in range(len(possible_labels)):
                    if len(summing_list) == 0:
                        np.append(summing_list, possible_labels[j])
                    elif len(summing_list) > 0 and summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                        summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + \
                                                                 possible_labels[j][0]
                    else:
                        np.append(summing_list, possible_labels[j])
                #summing_list = sorted(summing_list, key=operator.itemgetter(0))
                summing_list = summing_list [summing_list [:, 0].argsort()]
                data[i][1] = summing_list[len(summing_list) - 1][1]
                for k in range(3, len_line):
                    np.delete(data, [i,3])



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
                            np.append(data[data[i][2][k][0]],new_tmp)

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

