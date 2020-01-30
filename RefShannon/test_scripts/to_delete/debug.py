import os, pdb, sys

#st = 'comp %s\tnode %s\tm=%d\tn=%d\tnum_known_paths=%d\tsparsity=%d\n'
fld = 'dmp/chr1_removeDupP/'
files = os.listdir(fld)

stat = [] #list of [m, n, known_path, sparsity]
comp_stat = {} #key - comp, val - [num of nodes, avg kp/(m*n), avg sp/(m*n)]

i = 0; j=0; T=len(files)/100;
#pdb.set_trace()
for f in files:
    i+=1
    if i>=T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed'%j); sys.stdout.flush()

    if f[0:5]!='debug': continue    

    with open('%s/%s'%(fld, f), 'r') as f2:

        isFirst = True
        for line in f2:

            tokens = line.split()
            comp = tokens[0]
            name = tokens[1]
            m = int(tokens[2])
            n = int(tokens[3])
            p = int(tokens[4])
            s = int(tokens[5])

            stat.append([comp,name, m,n,p,s])

            if isFirst==True:
                comp_stat[comp] = [0, 0, 0]
                isFirst = False

            comp_stat[comp][0]+=1
            comp_stat[comp][1]+= float(p)/(m*n)
            comp_stat[comp][2]+= float(s)/(m*n)

        comp_stat[comp][1]=comp_stat[comp][1]/comp_stat[comp][0]
        comp_stat[comp][2]=comp_stat[comp][2]/comp_stat[comp][0]

pdb.set_trace()