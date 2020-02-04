from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

from __future__ import division
from cvxopt import matrix,spmatrix, solvers, spdiag, sparse
import copy,random, numpy
import pdb
import sys, time

from RefShannon.path_decompose_sparse import path_decompose
#from path_decompose_modi.path_decompose_sparse2 import path_decompose2 #change to: output = path_decompose2(
#from path_decompose_modi.path_decompose_sparse4 import path_decompose4 #change to: output = path_decompose4(
from RefShannon.path_decompose_modi.path_decompose_sparse5 import path_decompose5

#from memory_profiler import profile
#memory_fp = open('gen_sf_memory_profiler_non_single_nodes.log', 'w+')

modi_zero_edge_flow = False #Recommend: False
                            #When True: if n1-->v-->n2, n1-->v has 0 in-flow,  connect n1 to end;
                            #                           v-->n2 has 0 out-flow, connect start to n2;

import os
from time import sleep
from util import run_cmd
import numpy as np
#from dump_sf import describe_one_sample

sys.setrecursionlimit(100000)
solvers.options['show_progress'] = False
solvers.options['msg_lev'] = 'GLP_MSG_OFF'
use_norm = 'l1' #Options 'l1' or 'l2'
unit_normalization = False
restored_normalization = True
use_Y_paths = True
use_smoothing = False
use_GLPK = False
path_sparsity = 10

burden_factor = 100

run_penalized = 0

debug_mode = 0
overwrite_normalization = 0 # Set this to 1 if you want to overload the variable normalization to the copy_count of the true
#if overwrite_normalization=1: equivalent to using original Copy counts for thresholding
#pdb.set_trace()

dump = False #True

'''
usage:
python algorithm_SF.py comp -I intermediate_dir [-O algo_out_dir] [-tr_name_tag tr_name_tag] [-target target] [-F F_val] [--outputFasta]
'''
args = sys.argv
comp = args[1]
sample_name = args[args.index('-I')+1]

if '-O' in args:
    sample_name_out = args[args.index('-O')+1]
else:
    sample_name_out = parent_dir(sample_name)+'/algo_output/'

if '-tr_name_tag' in args:
    tr_name_tag = args[args.index('-tr_name_tag')+1]
else:
    tr_name_tag = 'refShannon'

if '-target' in args: #used for GTF output
    target = args[args.index('-target')+1]
else:
    target = '.'

if '-F' in args:
    F_val = float(args[args.index('-F')+1])
else:
    F_val = 0.0

outputGTF = True

outputFasta = True #default False
if '--outputFasta' in args:
    outputFasta = True

edges_file = sample_name+'/edges' + comp + '.txt'
nodes_file = sample_name+'/nodes' + comp + '.txt'
single_nodes_file = sample_name+'/single_nodes.txt'
KnownPathsFile = sample_name+'/paths' + comp + '.txt'
reconstr_file = sample_name_out+'/reconstructed_comp_' +str(comp) + '.fasta'
reconstr_Y_file = sample_name_out+'/reconstructed_comp_' +str(comp) + '.fasta' # 'algo_output/reconstructed' + '.fasta'

if dump == True:
    debug_file = sample_name_out+'/debug_comp_' +str(comp) + '.dmp'

if outputGTF == True:
    comp_gtf_file = sample_name_out+'/reconstructed_comp_' +str(comp) + '.gtf' # 'algo_output/reconstructed' + '.fasta'

use_file =1 #whether you want to read from files or not

if not os.path.exists(sample_name_out):
    run_cmd('mkdir -p ' + sample_name_out)

hash_to_node = {}      ## Map: hash value => corresponding node.
node_to_hash = {}      ## Map: node => node's hash value.
known_paths = []       ## A list of known paths, each path is a tuple of three elements
known_paths_str = []
paths_for_node = {}    ## Map: node => known_paths that the node is involved in.

#@profile(stream=memory_fp)
def single_nodes_to_fasta(): 
    ## The function outputs individual nodes without any edges as reconstructed transcripts.

    #pdb.set_trace()
    if outputFasta == True:
        reconstFile = open(reconstr_file, 'w')
    if outputGTF == True:
        gtf_file = open(comp_gtf_file, 'w') #'a'-->'w'

    i = 0
    for lines in open(single_nodes_file):
        if i>0:
            fields = lines.strip().split()
            #reconstFile.write('>Shannon_'+sname + '_single_'+str(i)+'\t Copycount:' + fields[2])
            transcript_id = tr_name_tag+'_single_node_tr_'+str(i)

            if outputFasta==True:
                reconstFile.write('>'+transcript_id+'\tweight=' + fields[2])
                reconstFile.write('\n'+fields[1] +'\n')
            #pdb.set_trace()
            if outputGTF==True:
                #pdb.set_trace()
                path_regions = '%d, %d;'%(int(fields[4]), int(fields[5])) #0-based, inclusive
                write_gtf(gtf_file, target, transcript_id, path_regions)
        i += 1

    if outputFasta == True:
        reconstFile.close()
    if outputGTF == True:
        gtf_file.close()
    #pdb.set_trace()


#@profile(stream=memory_fp)
def ParseKnownPathsFile(KnownPathsFile, graph): 
    ## This function builds the known_paths and paths_for_node data structures.
    f = open(KnownPathsFile, 'r')
    lines = f.readlines()
    i = 0
    for node1 in graph.nodes:
        paths_for_node[node1] = []
    for (i,line) in enumerate(lines):
        if i != 0:
            tokens = line.split()
            nodes_in_path = []
            #tmp_string = ""
            #prev_node = None
            for (j,hashcode) in enumerate(tokens):
                node = hash_to_node[hashcode]
                nodes_in_path.append(node); # pdb.set_trace()
                if paths_for_node.get(node) == None:
                    paths_for_node[node]=[i-1]
                else:
                    if len(paths_for_node[node]) < path_sparsity: #Append only if the no of known paths is smaller than path_sparsity
                        paths_for_node[node].append(i-1)
                #prev_node = node
            known_paths.append(nodes_in_path); #pdb.set_trace()
    f.close()


# must be called first
#@profile(stream=memory_fp)
def ParseNodeFile(NodeFile, graph):
    ## Builds node_to_hash and hash_to_node.
    ## NodeFile: The file with all the nodes.
    ## graph: the graph object we are using.
    f=open(NodeFile,'r')
    #lines = f.readlines()
    i = 0
    for line in f: #lines:
        if i != 0:
            tokens = line.split()
            try:
                t2=float(tokens[2])
            except ValueError:
                t2 = 0
            try:
                t3=float(tokens[3])
            except ValueError:
                t3 = 0

            #if tokens[1]=='*':
            #    pdb.set_trace()
            #    outputFasta = False

            if outputFasta == False:
                tokens[1] = ""

            if outputGTF==True:
                #pdb.set_trace()
                stt = int(tokens[4])
                stp = int(tokens[5])
                new_node = Node(tokens[1], t2,t3,tokens[0], stt, stp) #bases, cpcnt, norm, id, stt, stp
            else:
                new_node = Node(tokens[1], t2,t3,tokens[0]) #bases, cpcnt, norm, id, stt, stp

            hash_to_node[tokens[0]] = new_node
            node_to_hash[new_node] = tokens[0]
            graph.add_node(new_node)
        i += 1
    f.close()

#@profile(stream=memory_fp)
def ParseEdgeFile(EdgeFile, graph):
    ## Adds each edge to list of connections for both nodes involved.  *EHC
    ## EdgeFile: The file with all the edges.
    ## graph: the graph object we are using.
    f = open(EdgeFile, 'r')
    #lines = f.readlines()
    i = 0
    for line in f: #lines:
        if i != 0:
            tokens = line.split()
            start_node = hash_to_node[tokens[0]]
            end_node = hash_to_node[tokens[1]]
        
            start_node.out_edges.append([end_node, int(tokens[2]), float(tokens[3]), float(tokens[4])])
            end_node.in_edges.append([start_node, int(tokens[2]), float(tokens[3]), float(tokens[4])])
        i += 1
    f.close()

#@profile(stream=memory_fp)
def write_gtf(gtf_file, target, transcript_id, path_regions):

    path_regions_1 = [[int(itm.split(',')[0]), int(itm.split(',')[1])] for itm in path_regions.split(';') if itm != '' and '-1' not in itm] # list of [stt, stp]
    path_regions_2 = []
    
    #merge
    e_merge = []
    for e in path_regions_1:
        if e_merge == []:
            e_merge = e
        elif e_merge[-1]+1==e[0]:
            e_merge = [e_merge[0], e[1]]
        else:
            path_regions_2.append(e_merge)
            e_merge = e

    if e_merge != []:
        path_regions_2.append(e_merge)

    if path_regions_2==[]: return

    st =  '%s\t'%target
    st += 'refShannon\t'
    st += 'transcript\t'
    st += '%d\t'%(path_regions_2[0][0]+1) #1-based, inclusive
    st += '%d\t'%(path_regions_2[-1][1]+1)
    st += '.\t' #score
    st += '+\t' #strand
    st += '.\t' #frame
    st += 'transcript_id \"%s\";\n'%transcript_id #attribute

    if len(path_regions_2)>=1: #add == 1 because needed when use gffread to convert gtf to fa

        for i in range(len(path_regions_2)):
            st +=  '%s\t'%target
            st += 'refShannon\t'
            st += 'exon\t'
            st += '%d\t'%(path_regions_2[i][0]+1) #1-based, inclusive
            st += '%d\t'%(path_regions_2[i][1]+1)
            st += '.\t' #score
            st += '+\t' #strand
            st += '.\t' #frame
            st += 'transcript_id \"%s\"; exon_number \"%d\";\n'%(transcript_id, i+1) #attribute

    gtf_file.write(st)

    return

def intersect(a, b, c ):
    ## Returns the intersection of the three lists.  
     return list(set(a) & set(b) & set(c))
     
def intersect5(a, b, c, d, e):
    ## Returns the intersection of the five lists.  
     return list(set(a) & set(b) & set(c) & set(d) & set(e))
     
class Edge(object):  ## Edge object (used for building copycount filter matrices).
    def __init__(self, start_node, end_node, overlap_weight, weight, L):
        self.start = start_node  ## The starting node in the edge.
        self.end = end_node  ## The ending node in the edge.
        self.overlap_weight = overlap_weight  ## How much the connected nodes overlap.
        self.weight = weight  ## The copycount for the edge.
        self.L = L  ## The normalization used for error term in copycount filtering. 

class Node(object):  ## Node object (used universally)
    def __init__(self, node_string, node_weight, L,name, stt=-1, stp=-1):
        self.string = node_string  ## The sequence of bases that the node represents.
        self.in_edges = []  ## A list of edges in which this node object is the end node.
        self.out_edges = []  ## A list of edges in which this node object is the start node.
        self.name = name  ## Hash value for the node.
        self.weight = (node_weight)  ## The copycount for the node.
        self.L = (L)  ## The normalization used for the error term in copycount filtering.  
        self.DNA_start_pos = stt  ## Start position in the reference DNA at which this node has bases from. 0-based, inclusive
        self.DNA_end_pos = stp  ## Start position in the reference DNA at which this node has bases from. 0-based, inclusive
    def set_string(self, node_string):  ## This is the sequence on the node.
        self.string = node_string
    def add_in_edge(self, node, overlap_weight, weight, L):
        self.in_edges.append([node, overlap_weight, weight, L])
    def add_out_edge(self, node, overlap_weight, weight, L):
        self.out_edges.append([node, overlap_weight, weight, L])

class Graph(object):    ## Graph object (used universally)
    def __init__(self):
        self.start = None    ## Node with no in-edges
        self.end = None      ## Node with no out-edges
        self.nodes = []
        self.tobereduced = []  ## list of nodes with more than one in-edge
        self.paths = []
        #for matrix
        self.edges = []    ## This is used for building the matrix with the node's in/out edge information.
        #for edge filter
        self.edges2 = []
        self.edge_weights = []
        self.node_weights = []
        self.normalization = []
        self.penalization = []
        #for unique solution determination
        self.no_unique_solution = False      ## Is there a unique sparsest flow for each decomposed node in the graph?
        self.paths=[]       ## known paths in graph
        self.paths_Y = []   
        self.og_nodes = {}  ## original nodes in graph before sparse flow is run
        self.constituent_nodes = {}  ## dictionary with nodes as keys and nodes that were condensed to form node as values


    def add_node(self, node):
        self.nodes.append(node)
    def remove_node(self, node):
        self.nodes.remove(node)
    def add_edge(self, start_node, end_node, overlap_weight, weight, L):
        start_node.out_edges.append([end_node, overlap_weight, weight, L])
        end_node.in_edges.append([start_node, overlap_weight, weight, L])

    def findStartAndEnd2(self):
        ## Adds a dummy node labeled start_node that has an out edge to all original nodes with in degree 0.
        ## Adds a dummy node labeled end_node that has an in edge to all original nodes with out degree 0.  
        start_node = Node("Start_", 0, 0,'S')
        end_node = Node("_End", 0, 0,'E')

        for node in self.nodes:
            if len(node.in_edges) == 0:
                node.in_edges.append([start_node, 0, node.weight, 0])
                start_node.out_edges.append([node, 0, node.weight, 0])
                start_node.weight += float(node.weight)

            if len(node.out_edges) == 0:
                node.out_edges.append([end_node, 0, node.weight, 0])
                end_node.in_edges.append([node, 0, node.weight, 0])
                end_node.weight += float(node.weight)

        self.nodes.append(start_node)
        self.nodes.append(end_node)
        self.start = start_node
        self.end = end_node


    def printNodes(self):
        ##  Prints out each node with all in-edges and out-edges on the same line.
        print('Nodes:\n' , [[e.name, e.weight, e.L ] for e in self.nodes])
        print('\n')
        for each in self.nodes:
            if len(each.out_edges) == 0:
                list1 = [[each.in_edges[i][0].name,  each.in_edges[i][1], each.in_edges[i][2], each.in_edges[i][3]] for i in range(0, len(each.in_edges))]
                print(each.name,"   out edges: None", "   in edges:",  list1)
            if len(each.in_edges) == 0:
                list1 = [[each.out_edges[i][0].name,  each.out_edges[i][1], each.out_edges[i][2], each.out_edges[i][3]] for i in range(0, len(each.out_edges))]
                print(each.name,"   out edges:", list1, "   in edges: None")
            if len(each.out_edges) != 0 and len(each.in_edges) != 0:
                list_out = [[each.out_edges[i][0].name,  each.out_edges[i][1], each.out_edges[i][2], each.out_edges[i][3]] for i in range(0, len(each.out_edges))]
                list_in = [[each.in_edges[i][0].name,  each.in_edges[i][1], each.in_edges[i][2], each.in_edges[i][3]] for i in range(0, len(each.in_edges))]
                print(each.name,"   out edges:", list_out, "in edges:",list_in)


    def printNodesSmall(self):
        ##  Prints out each node with with it's out edges on the same line.
        for each in self.nodes:
            if len(each.out_edges) != 0:
                list_out = [[each.out_edges[i][0].name] for i in range(0, len(each.out_edges))]
                print(each.name,"   out edges:", list_out)

    def search(self):
        '''Searches for all nodes in the graph that have more than one in edge OR more than one out edge and adds them to
        to the list of nodes to be reduced by path_decompose.  
        If use_Y_paths is on, only sends nodes to path_decompose if the node has one in edge and more than one out edge
        Sorts nodes in topological order.
        '''
        del self.tobereduced[:]
        for node in self.nodes:
            if (len(node.in_edges)!=0)  and (len(node.out_edges)!=0) and (node is not self.start) and (node is not self.end):
                if len(node.in_edges) >1 or len(node.out_edges)>1: #SEARCH FOR any Y nodes
                    if use_Y_paths and len(node.in_edges)<=1:  #Search for left-Y NODES and X-nodes
                        continue
                    self.tobereduced.append(node)
        self.tobereduced.sort(key=lambda x: int(x.name.split("_")[0]), reverse=False)

    
    #@profile(stream=memory_fp)
    def algorithm2(self):
        # Runs sparse flow algorithm on graph to simplify graph such that all nodes have an in-degree of 1.
        done = False
        cycle_limit = 1
        algo_iteration= 0
        node_iter = 0 

        for each in self.nodes:
            self.constituent_nodes[each] = [each]
        ## If wau == True, you need to worry about condensing.
        while not done:
            #new_node_list = copy.copy(self.nodes) # Algorithm will construct new list of nodes for each iteration through all nodes
            
            for node in self.nodes:
                if node == self.start or node == self.end:
                    continue
                if use_Y_paths or algo_iteration == 0:
                    if len(node.in_edges) <= 1:
                        continue
                else:
                    if len(node.in_edges)<=1 and len(node.out_edges)<=1:
                        continue
                node_iter += 1; 
                #sys.stdout.write('\r')
                #sys.stdout.write(str(time.asctime())+': comp: ' + str(comp) + ', algo_iter: ' + str(algo_iteration)  + ', node_iter: ' +  str(node_iter) + ', node_name: ' +str(node.name)+ ', m: ' + str(len(node.in_edges)) + ', n: ' + str(len(node.out_edges)) + ', Paths: ' + str(len(paths_for_node.get(node,[]))))
                #sys.stdout.flush(); 
                if 1:
                    if 1:
                        new_nodes = []  ## list of new nodes produced from decomposition of the current node.
                        inedges = []  ## A vector that contains the connected node of each in-edge.
                        outedges = []  ## A vector that contains the connected node of each out-edge.
                        inedge_vector = []  ## A vector that contains the copycounts of each in-edge.
                        outedge_vector = []  ## A vector that contains the copycounts of each in-edge.
                        inedge_cc = []  ##  A vector that should contain the copycounts of each in node, but currently contains the overlap of the sequence.
                        outedge_cc = []  ##  A vector that should contain the copycounts of each out node, but currently contains the overlap of the sequence.


                        incoming_edge_attributes = {}  ## A dictionay that contains the overlap and normalization information for each in edge.
                        outgoing_edge_attributes = {}  ## A dictionay that contains the overlap and normalization information for each out edge.
                        if len(node.in_edges) == 0:
                            #print('Hanging Node!'); 
                            node.in_edges.append([self.start, 0, node.weight, 0])
                            self.start.out_edges.append([node, 0, node.weight, 0])
                            self.start.weight += float(node.weight)

                        if len(node.out_edges) == 0:
                            #print('Hanging Node!'); 
                            node.out_edges.append([self.end, 0, node.weight, 0])
                            self.end.in_edges.append([node, 0, node.weight, 0])
                            self.end.weight += float(node.weight)
 
                        for (i,in_edge) in enumerate(node.in_edges):
                            inedges.append(in_edge[0])
                            inedge_vector.append(float(in_edge[2]))
                            inedge_cc.append(float(in_edge[3]))
                            incoming_edge_attributes[in_edge[0]] = [in_edge[1], in_edge[3]]

                        for (j,out_edge) in enumerate(node.out_edges):
                            outedges.append(out_edge[0])
                            outedge_vector.append(float(out_edge[2]))
                            outedge_cc.append(float(in_edge[3]))
                            outgoing_edge_attributes[out_edge[0]] = [out_edge[1], out_edge[3]]

                        P = matrix(0.,(len(node.in_edges), len(node.out_edges)))
            
                        #  This section of code determines which known paths will be considered when decomposing this node
                        # pdb.set_trace()
                        #tmp_paths = []
                        #'''
                        '''in_node_cnt = {} #key - in_node_org, val - cnt (e.g. two innodes can come from same org node)
                        out_node_cnt = {}
                        for (m, in_node) in enumerate(inedges):
                            in_node_org = self.constituent_nodes[in_node][0]
                            if in_node_org in in_node_cnt:
                                in_node_cnt[in_node_org] += 1; #pdb.set_trace()
                            else:
                                in_node_cnt[in_node_org] = 1

                        for (n, out_node) in enumerate(outedges):
                            out_node_org = self.constituent_nodes[out_node][0]
                            if out_node_org in out_node_cnt:
                                out_node_cnt[out_node_org] += 1; #pdb.set_trace()
                            else:
                                out_node_cnt[out_node_org] = 1'''

                        I0 = set(); J0 = set()

                        for (m, in_node) in enumerate(inedges):
                            for (n, out_node) in enumerate(outedges):
                                if '_' in in_node.name or '_' in out_node.name:
                                    #pdb.set_trace()
                                    continue #less P=1
                                in_node_org = self.constituent_nodes[in_node][0]
                                c_node_org = self.constituent_nodes[node][0]
                                out_node_org = self.constituent_nodes[out_node][0]
                                if [in_node_org, c_node_org, out_node_org] in known_paths:
                                    # and \
                                    #in_node_cnt[in_node_org]<=1 and \
                                    #out_node_cnt[out_node_org]<=1:
                                    P[m,n]=1; #pdb.set_trace()

                                    I0.add(m); J0.add(n)

                        if 0: #algo_iteration==0:
                            if (len(I0)==m and len(J0)==n) or (len(I0)==m-1 and len(J0)==n-1):
                                #pdb.set_trace() #LS approach
                                output = path_decompose5(inedge_vector, outedge_vector, inedge_cc, outedge_cc, overwrite_normalization, P,use_GLPK, path_sparsity)
                            else:
                                #pdb.set_trace() #next node
                                continue 
                        else:
                            #'''
                            #            tmp_paths.append([m,n])
                            #if len(inedges)>1 and len(outedges)>1 and tmp_paths!=[]:
                            #    print('m:%d'%len(inedges))
                            #    print('n:%d'%len(outedges))
                            #    print('%s\n'%(str(tmp_paths)))
                            #    #pdb.set_trace()                 
                            #  This line decomposes the node 
                            #if sum([1 for i in inedge_vector if i<0])>0 \
                            #   or sum([1 for j in outedge_vector if j<0])>0:
                            #    pdb.set_trace()               
                            #output = path_decompose(inedge_vector, outedge_vector, inedge_cc, outedge_cc, overwrite_normalization, P,use_GLPK, path_sparsity)
                            output = path_decompose(inedge_vector, outedge_vector, inedge_cc, outedge_cc, overwrite_normalization, P,use_GLPK, path_sparsity, F_val)
                        temp_matrix = output[0]
                        m = len(inedge_vector)
                        n = len(outedge_vector)
                        #pdb.set_trace()

                        '''if modi_zero_edge_flow==True:
                            in_node_flow = numpy.sum(temp_matrix, 1)
                            I_no_flow = [inedges[i] for i in range(m) if in_node_flow[i]==0] #list of in-nodes {r} s.t. no flow decomposed through r-->v
                            out_node_flow = numpy.sum(temp_matrix, 0)
                            J_no_flow = [outedges[j] for j in range(n) if out_node_flow[j]==0] '''#list of out-nodes {t} s.t. no flow decomposed through v-->t
                            #if I_no_flow != [] or J_no_flow != []:
                            #    pdb.set_trace()
                        
                        nodes_to_eliminate = [node]
                        
                        # This section of the code builds the new nodes formed during decomposition, and implicitly condenses the 1x1 nodes
                        for i in range(0, m):
                            for j in range(0, n):
                                curr_edge_cc = temp_matrix[i][j]
                                if curr_edge_cc != 0:
                                    out_attr = outgoing_edge_attributes[outedges[j]]  
                                    in_attr = incoming_edge_attributes[inedges[i]]
                                    if outputGTF == True:
                                        new_node = Node(node.string, curr_edge_cc, node.L,node.name+"_["+str(i)+","+str(j)+"]", node.DNA_start_pos, node.DNA_end_pos)
                                    else:
                                        new_node = Node(node.string, curr_edge_cc, node.L,node.name+"_["+str(i)+","+str(j)+"]")
                                    new_node.add_in_edge(inedges[i], in_attr[0], curr_edge_cc, in_attr[1])
                                    inedges[i].add_out_edge(new_node,in_attr[0], curr_edge_cc, in_attr[1])
                                    new_node.add_out_edge(outedges[j], out_attr[0], curr_edge_cc, out_attr[1])
                                    outedges[j].add_in_edge(new_node, out_attr[0], curr_edge_cc, out_attr[1])
                                    self.nodes.append(new_node)
                                    self.constituent_nodes[new_node] = self.constituent_nodes[node]

                        '''
                        with open('dmp/dmp_P_info.txt', 'a') as f_dmp_P_info:
                            st = 'comp %s\tnode %s\tm=%d\tn=%d\n'%(comp, node.name, len(inedges), len(outedges))
                            st += describe_one_sample(temp_matrix, len(inedges), len(outedges), P, inedge_vector, outedge_vector)
                            f_dmp_P_info.write(st)
                        '''
                        if dump == True:
                            with open(debug_file, 'a') as dbg_f:
                                #pdb.set_trace()
                                nonzero_P = np.count_nonzero(P)
                                nonzero_res = np.count_nonzero(temp_matrix)
                                #st = 'comp %s\tnode %s\tm=%d\tn=%d\tnum_known_paths=%d\tsparsity=%d\n'%(comp, node.name, \
                                #                                                                        len(inedges), len(outedges), \
                                #                                                                        nonzero_P, nonzero_res)
                                st = '%s\t%s\t%d\t%d\t%d\t%d\n'%(comp, node.name, \
                                                                 len(inedges), len(outedges), \
                                                                 nonzero_P, nonzero_res)
                                dbg_f.write(st)



                        ## For each node that was condensed into a new node, delete all it's connections.
                        for edge in node.in_edges:
                            in_node_temp = edge[0]
                            for oedge in in_node_temp.out_edges:
                                if oedge[0] is node:
                                #if oedge[0].string == node.string:
                                    in_node_temp.out_edges.remove(oedge)

                            '''if modi_zero_edge_flow==True \
                               and in_node_temp in I_no_flow and len(in_node_temp.out_edges)==0 \
                               and in_node_temp != self.start:
                                #pdb.set_trace()
                                in_node_temp.out_edges.append([self.end, 0, in_node_temp.weight, 0])
                                self.end.in_edges.append([in_node_temp, 0, in_node_temp.weight, 0])
                                self.end.weight += float(in_node_temp.weight)'''

                        for edge in node.out_edges:
                            out_node_temp = edge[0]
                            for iedge in out_node_temp.in_edges:
                                if iedge[0] is node:
                                    out_node_temp.in_edges.remove(iedge)

                            '''if modi_zero_edge_flow==True \
                               and out_node_temp in J_no_flow and len(out_node_temp.in_edges)==0 \
                               and out_node_temp != self.end:
                                #pdb.set_trace()
                                out_node_temp.in_edges.append([self.start, 0, out_node_temp.weight, 0])
                                self.start.out_edges.append([out_node_temp, 0, out_node_temp.weight, 0])
                                self.start.weight += float(out_node_temp.weight)   '''

                        if node not in self.nodes:
                            'alert'
                        else:    
                            self.nodes.remove(node)
                        #if node not in self.constituent_nodes:
                        #    'alert2'
                        #else:
                        #    #pdb.set_trace()
                        #    del self.constituent_nodes[node]; #pdb.set_trace()
            #self.nodes = new_node_list # update node list after each iteration through all nodes
            self.search() # checks to see if any more nodes need to be reduced

            if 0: #algo_iteration==0:

                #pdb.set_trace()

                # This is to ensure nodes are run through topologically
                self.nodes.remove(self.start)
                self.nodes.remove(self.end)
                self.nodes.sort(key=lambda x: int(x.name.split("_")[0]), reverse=False)
                self.nodes.append(self.end)
                self.nodes.insert(0, self.start)

            else: #2nd iteration, ready to stop

                if len(self.tobereduced) == 0:
                    done = True
                else:
                    # This is to ensure nodes are run through topologically
                    self.nodes.remove(self.start)
                    self.nodes.remove(self.end)
                    self.nodes.sort(key=lambda x: int(x.name.split("_")[0]), reverse=False)
                    self.nodes.append(self.end)
                    self.nodes.insert(0, self.start)
                
            algo_iteration += 1
        #sys.stdout.write('\n')               

    def read_paths_recursive_modi(self, node, str_till_now, nodes_till_now, regions_till_now, overlap,prev_weight):
        '''Reads all paths in graph recursively
        node:  Current node
        str_till_now:  The string seen before this node.
        overlap:  THe amount of bases of overlap between the last node in the path and the current node.
        prev_weight:  The wieght of thw last node in the path.  
        '''
        # pdb.set_trace()
        if outputFasta==False:
            curr_str = str_till_now + '*' # *** (no START_)
        else:
            curr_str= str_till_now + node.string[overlap:]
        node_name = node.name.split('_')[0]
        curr_nodes = nodes_till_now + '->'+ node_name
        curr_regions = regions_till_now + '%d, %d;'%(node.DNA_start_pos, node.DNA_end_pos)
        if len(node.out_edges) == 0: ## This assumes all paths end at the _END node.
            #if curr_str[-4:] != '_End':
            if node.string != '_End':
                #Return without appending this path
                return
            else:
                if outputFasta==True: 
                    curr_str = curr_str[:-4] #START_...(no _End) or *** 
                self.paths_Y.append([curr_str,prev_weight,curr_nodes, curr_regions])
                return
        else:
            prev_weight = node.weight
            for (i,each) in enumerate(node.out_edges):
                new_node=each[0]
                overlap = int(each[1])
                #pdb.set_trace()
                self.read_paths_recursive_modi(new_node,curr_str,curr_nodes,curr_regions, overlap,prev_weight)


    def read_Y_paths_modi(self):
        ''' Uses read_paths_recursive to find all paths if the graph only has Y nodes 
        (a Y node is a node with at most 1 in edge AND 0 or more out edges).
        '''
        if outputFasta == True:
            pathfile = open(reconstr_Y_file, 'w')
        if outputGTF == True:
            gtf_file = open(comp_gtf_file, 'w')

        self.search()
        if len(self.tobereduced) != 0:
            print('CAUTION:There are still some unresolved nodes')
        self.read_paths_recursive_modi(self.start,'','','', 0,0)
        #pdb.set_trace()
        for (i,path_info) in enumerate(self.paths_Y): #self.paths_Y.append([curr_str,prev_weight,curr_nodes, curr_regions])
            if outputFasta==True:
                path_str = path_info[0][6:]#skip START_
            else:
                path_str = path_info[0]
            path_wt = path_info[1]
            path_nodes = path_info[2]
            path_regions = path_info[3]
            
            if len(path_str):# *** (if outputFasta==False) or [ATCG]+ (if outputFasta==True)
                transcript_id = tr_name_tag+'_comp_'+comp+'_tr_'+str(i)
                if outputFasta==True:
                    pathfile.write('>'+transcript_id+"\tweight="+str(path_wt)+'\tnodes='+path_nodes)
                    pathfile.write("\n"+path_str+"\n") #with weights
                if outputGTF==True:
                    write_gtf(gtf_file, target, transcript_id, path_regions)

        if outputFasta == True:
            pathfile.close()
        if outputGTF == True:
            gtf_file.close()

#main

#pdb.set_trace()
if comp == '-1':
    single_nodes_to_fasta()
    sys.exit(0)

graph2 = Graph()
ParseNodeFile(nodes_file, graph2)
ParseEdgeFile(edges_file, graph2)
ParseKnownPathsFile(KnownPathsFile, graph2)

graph2.findStartAndEnd2()

if len(graph2.nodes) <= 3:
    graph2.read_Y_paths_modi()
    sys.exit(0)
    
#t_start = time.time()
graph2.algorithm2()
#t_elapsed = (time.time() - t_start)

graph2.read_Y_paths_modi()










                                    
                                    
                                    
                                    

