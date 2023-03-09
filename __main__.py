import argparse
import logging
import sys
from igraph import *
from labprop.LabelPropagation import lp1
import time
from multiprocessing import Pool
#from numba import njit
import numpy as np
import networkx as nx


# Sample command
# -------------------------------------------------------------------
# python ReadGraph_SGA.py     --graph /path/to/graph_file.asqg
#                            --binned /path/to/binning_result.csv
#                            --output /path/to/output_folder
#                            --max_iteration maximum number of iterations
#                            --read_type 1 paired-end 2 single-end
#                            --assembler 1 sga 2 minimap2
# -------------------------------------------------------------------

if __name__ == "__main__":

    # Setup logger
    # -----------------------

    logger = logging.getLogger('ClassGraph')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    logger.addHandler(consoleHeader)

    # Setup argument parser
    # -----------------------

    ap = argparse.ArgumentParser()

    ap.add_argument("--graph", required=True, help="path to the reads graph file")
    ap.add_argument("--binned", required=True,
                    help="path to the file with the initial binning output")
    ap.add_argument("--output", required=True, help="path to the output folder")
    ap.add_argument("--prefix", required=True, help="prefix for the output file")
    ap.add_argument("--max_iteration", required=False, type=int, default=20,
                    help="maximum number of iterations for label propagation algorithm. [default: 20]")
    #ap.add_argument("--lp_version", required=False, type=int, default=1,
    #                help="Type 1 if you want to propagate with lp-v1, type 2 if you prefer to use lp-v2. [default 1]")
    ap.add_argument("--read_type", required=False, type=int, default=1,
                    help="Type 1 if the reads are paired, type 2 if reads are single-end[default 1]")
    ap.add_argument("--assembler", required=False, type=int, default=1,
                    help="Type 1 if the assembler is sga, type 2 if is minimap2")
    #ap.add_argument("--beta", required=False, type=float, default=0,
    #                help="Choose the value of Beta between 0 and 1 ( 1 excluded )")

    args = vars(ap.parse_args())

    sgafile = args["graph"]
    # asqg
    kraken2_file = args["binned"]
    # kraken
    output_path = args["output"]
    prefix = args["prefix"]
    max_iteration = args["max_iteration"]
    #labprop_v = args["lp_version"]
    read_type = args["read_type"]
    assembler = args["assembler"]
    #beta = args["beta"]


    # Setup output path for log file
    # ---------------------------------------------------

    fileHandler = logging.FileHandler(output_path + "/" + prefix + "classgraph.log")
    fileHandler.setLevel(logging.INFO)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info("ReadGraph makes use of the assembly graph produced by the assembler")
    logger.info("Overlap graph file: " + sgafile)
    logger.info("Existing binning output file: " + kraken2_file)
    logger.info("Final binning output file: " + output_path)
    logger.info("Maximum number of iterations: " + str(max_iteration))
    #logger.info("Label propagation v" + str(labprop_v))
    # logger.info("Beta = " + str(beta))
    # if beta > 1:
    #     logger.error("The value of beta must be in the interval [0,1]")
    #     logger.info("Exiting ClassGraph... Bye...!")
    #     sys.exit(1)
    logger.info("ReadGraph started")

    # Get the number of bins from the initial binning result
    # --------------------------------------------------------

    logger.info("Combining the graph file with the classification one, in order to obtain a labelled graph")

    krakenlabels = []
    links = []

    read_dict = {}
    inv_read_dict = {}

    reads_info = []

    links_dict = {}
    index = 0
    bool_idx1 = False
    bool_idx2 = False
    max_seq_len = 0

    def add_to_dicts(read_id):
        read_dict[index] = read_id
        inv_read_dict[read_id] = index

    def compute_normalize_overlap_len(assembler, line, max_seq_len):
        if assembler == 1:
            overlapstart = int(line.split()[3])
            overlapend = int(line.split()[4])
            normalizedoverlaplength = float((overlapend - overlapstart + 1) / max_seq_len)
        else:
            overlaplen = int(line.split()[11])
            normalizedoverlaplength = float(overlaplen / max_seq_len)
        return normalizedoverlaplength

    with open(sgafile, 'r') as graph_file, open(kraken2_file, 'r') as classification_file:
        # Build the following dictionary: {0: 'read_id1', 1: 'read_id2, ...}
        for line in classification_file:
            read_id = line.split()[0]
            tax_id = line.split()[1]
            read_dict[index] = read_id
            reads_info.append(tax_id)
            index += 1

        index -= 1

        inv_read_dict = {v: k for k, v in read_dict.items()}

        for line in graph_file:
            if line.startswith("ED"):
                seq1_len = int(line.split()[5])
                seq2_len = int(line.split()[8])
                max_seq_len = max(max_seq_len, seq1_len, seq2_len)

        graph_file.seek(0)

        for line in graph_file:
            if line.startswith("VT"):
                # If we have paired end read
                if read_type == 1:
                    read_id = line.split()[1][:-2]
                    if read_id not in inv_read_dict:
                        reads_info.append(0)
                        index += 1
                        add_to_dicts(read_id)
                else:
                    read_id = line.split()[1]
                    if read_id not in inv_read_dict:
                        reads_info.append(0)
                        index += 1
                        add_to_dicts(read_id)
            if line.startswith("ED"):
                if read_type == 1:
                    read1_id = line.split()[1][:-2]
                    read2_id = line.split()[2][:-2]
                else:
                    read1_id = line.split()[1]
                    read2_id = line.split()[2]
                if read1_id not in inv_read_dict:
                    reads_info.append(0)
                    index += 1
                    index_r1 = index
                    bool_idx1 = True
                    add_to_dicts(read1_id)
                if read2_id not in inv_read_dict:
                    reads_info.append(0)
                    index += 1
                    index_r2 = index
                    bool_idx2 = True
                    add_to_dicts(read2_id)
                if bool_idx1 == True:
                    bool_idx1 = False
                    in_dict_1 = int(index_r1)
                else:
                    in_dict_1 = int(inv_read_dict[read1_id])
                if bool_idx2 == True:
                    bool_idx2 = False
                    in_dict_2 = int(index_r2)
                else:
                    in_dict_2 = int(inv_read_dict[read2_id])
                normalizedoverlaplength = compute_normalize_overlap_len(assembler, line, max_seq_len)
                to_dict_key = str(in_dict_1) + " " + str(in_dict_2)
                if assembler == 1:
                    if to_dict_key in links_dict:
                        prev_val = links_dict[to_dict_key]
                        new_val = max(prev_val, normalizedoverlaplength)
                        links_dict[to_dict_key] = new_val
                    else:
                        links_dict[to_dict_key] = normalizedoverlaplength
                else:
                    links_dict[to_dict_key] = normalizedoverlaplength

    for key, value in links_dict.items():
        tmp = [int(key.split()[0]), int(key.split()[1]), value]
        links.append(tmp)

    # Construct the assembly graph
    # -------------------------------
    #
    try:

        start_time2 = time.process_time()
        # Create the graph
        reads_graph = Graph()

        # Create list of edges
        edge_list = []

        weights_list = []

        # Add vertices
        reads_graph.add_vertices(index + 1)

        #Name vertices
        for i in range(len(reads_graph.vs)):
            reads_graph.vs[i]["id"] = i

        # Iterate links
        for i in range(len(links)):
            # Remove self loops
            if links[i][0] != links[i][1]:
                edge_list.append((int(links[i][0]), int(links[i][1])))
                weights_list.append((float(links[i][2])))

        reads_graph.add_edges(edge_list)

        reads_graph.es["weight"] = weights_list

        print(f"time2 ={time.process_time()-start_time2}")

    except:
        logger.error("Please make sure that the correct path to the assembly graph file is provided.")
        logger.info("Exiting ClassGraph... Bye...!")
        sys.exit(1)

    logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

    logger.info("Preparing data for label propagation")

    # In the graph are not inserted edges that connect already labbeled nodes

    data = []
    LabbeledVertices = []
    reads_info = np.array(reads_info, dtype=np.uint32)

    # Label Delition
    to_elim = []
    for read in range(index +1):
        neighbours = reads_graph.neighbors(read)
        n = {}
        lab_reads = 0        
        for v in neighbours: 
            if reads_info[v] != 0:
                n[reads_info[v]] = 1 + n.get(reads_info[v], 0)
                lab_reads+=1

        coherent_lab =reads_info[read]
        if coherent_lab < lab_reads/2:
            to_elim.append(read)
            
    for read in to_elim: 
        reads_info[read] = 0


    for read in range(index + 1):

        neighbours = reads_graph.neighbors(read)
        neighs = []

        for neighbour in neighbours:
            if int(reads_info[neighbour]) == 0:
                #e_weight = reads_graph[read][neighbour]["weight"]
                e_weight = reads_graph.es[reads_graph.get_eid(read, neighbour)]["weight"]
                n = np.array([neighbour, e_weight])
                neighs.append(n)

        if len(neighs) == 0:
            data.append(np.array([]))
        else:
            neighs_numpy = np.stack(neighs, axis=0)
            data.append(neighs_numpy)
        

    try:

        logger.info("Starting label propagation")
        lp1(max_iteration, data, reads_info)
       
    except:
        logger.error("Please make sure that you inserted the correct parameter for the lp version (either 1 or 2)")
        logger.info("Exiting ClassGraph.. Bye...!")
        sys.exit(1)

    output_file = output_path + prefix + 'CG.res'

    with open(output_file, mode='w') as out_file:
        for i in range(len(data)):
            out_file.write(read_dict[i] + "\t" + str(reads_info[i]) + "\n")

    logger.info("You will find the output here: " + output_file )
    logger.info("***************Label propagation termined**************")
