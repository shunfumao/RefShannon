import pdb, heapq, copy
#from path_decompose_sparse import *

#'''
def local_dot(a,b):
    ## This function takes the dot product between two vectors.  
    n = len(a)
    sum=0
    for i in range(n):
        sum += a[i]*b[i]
    return sum


def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]
#'''

def path_decompose4(a,b,a_true,b_true,overwrite_norm,P,use_GLPK, sparsity= False):
    '''This function takes in a node's information in the attempt to decompose it into the lowest number of paths that
    accounts for the flow constraints.  
    Thhe algorithm uses many trials of a randomizaed optimization model and takes the best result seen.  
    **Do not use over_write_norm becuase a_true and b_true surrently have normalization instead of copycount. 
    a: a vector of the copytcounts of the in-edges for the node.
    b: a vector of the copycounts of the out-edges for the node.  
    a_true: a vector that should have the copycount values for the in-edges but does not.  (Don't Use)
    b_true: a vector that should have the copycount values for the out-edges but does not.  (Don't Use)
    decides whether or not to  use a_true and b_true.
    This is a matrix of in-edges versus out-edges that has a 0 if there is a known path and a 1 otherwise.  
    '''

    #mb_check = 1 #if this parameter is set to 1, if the data can be set using MB, then it
    from cvxopt import matrix, solvers, spmatrix, printing
    import numpy, copy
    from numpy import linalg as LA
    from numpy import array
    from numpy import zeros
    import operator
    solvers.options['msg_lev'] = 'GLP_MSG_OFF'
    solvers.options['show_progress'] = False
    m = len(a)  ## a is a vector of the current in edge copycount values.
    n = len(b)  ## b is a vector of the current out edge copycount values.
    sa = sum(a); sb = sum(b)

    ## Trivial case 
    if m==0 or n==0:
        return [[],0]
    if m==1:
        answer = array(matrix(b,(m,n)))
        return [answer,0]
    elif n==1:
        answer = array(matrix(a,(m,n)))
        return [answer,0]

    if sa<=0 or sb<=0:
        answer = array(matrix(0,(m,n)))
        return [answer,0]        

    ## Make all in flow values non-zero.
    '''for i in range(m):
        a[i]=max(a[i],1e-10)
    for j in range(n):
        b[j]=max(b[j],1e-10)'''

    ## If the flow in does not equal the flow out, make them equal.
    if sa>sb:
        const = (sa-sb)
        b = [k+const*k/sb for k in b ]
    else:
        const = (sb-sa)
        a = [k+const*k/sa for k in a ]    

    '''
    A = matrix(0.,(m+n-1,m*n))  ## A is used to enforce flow constraint on decomposition.
    p=matrix(0.,(m*n,1))  ## This vector tells whether there is a known path for the pair of nodes or not.
    
    for i in range(m):
        for j in range(n):
            A[i,j*m+i] = 1.;  #check indexing
            p[j*m+i] = 1-P[i,j];#p denotes if known path is *ABSENT*
    
    for j in range(n-1):  #Must be range(n) to get full mattrix
        for i in range(m):
            A[m+j,j*m+i] = 1.;
    z = a +b

    rhs = matrix(z)
    rhs = rhs[0:n+m-1,0]  ## this vector is used to enforce flow constraints 

    G = spmatrix(-1.,range(m*n),range(m*n))
    h = matrix(0.,(m*n,1))

    scale_factor = max(max(rhs),1e-100) * 0.01;
    '''

    '''########## 1st trial

    p = matrix(0., (m*n,1)) ## f_ij corresponds to p[jm+i]
    I0 = set() #i in I0 <--> exists j s.t. P[i,j]==1
    J0 = set() #j in J0 <--> exists i s.t. P[i,j]==1
    for i in range(m):
        for j in range(n):
            p[j*m+i] = 1-P[i,j] # p[jm+i]==1: known path for [i,j] is absent
            if P[i,j]==1:
                I0.add(i); J0.add(j)

    mask = matrix(0., (m*n,1)) ## f_ij==0 corresponds to mask[jm+i]==1
    mask_indice = [] #list of (i,j) where mask_ij (mask[j*m+i]) is 1 (f_ij 0)
    for i in range(m):
        for j in range(n):
            if P[i,j]==0 and (i in I0 or j in J0):
                mask[j*m+i]=1
                mask_indice.append((i,j))

    rhs = []
    #constraint from inedges
    I_w = []; J_w = [] #I/J indice constrained by edge weights w
    for i in range(m):
        prod = 1
        for j in range(n):
            prod = prod * mask[j*m+i]
        if prod==0:
            rhs.append(a[i])
            I_w.append(i)
    #constraint from outedges
    for j in range(n):
        prod = 1
        for i in range(m):
            prod = prod * mask[j*m+i]
        if prod==0:
            rhs.append(b[j])
            J_w.append(j)
    if len(I_w)==m and len(J_w)==n:
        del rhs[-1]; del J_w[-1] #avoid singular KKT matrix condition
    #constraint from mask
    for i in range(len(mask_indice)):
        rhs.append(0)

    rhs = matrix(rhs)

    #pdb.set_trace()

    A = matrix(0.0, (len(I_w)+len(J_w)+len(mask_indice), m*n)) # constraints * # of f_ij
    #A = matrix(0.0, (len(I_w)+len(J_w), m*n)) # constraints * # of f_ij
    row_id = 0
    for i in I_w:
        for j in range(n):
            A[row_id, j*m+i] = 1.0
        row_id += 1

    for j in J_w:
        for i in range(m):
            A[row_id, j*m+i] = 1.0
        row_id += 1

    for i,j in mask_indice:
        A[row_id, j*m+i] = 1.0
        row_id += 1

    scale_factor = max(max(rhs),1e-100) * 0.01;

    #pdb.set_trace()

    '''##########

    ########## 2nd trial

    p = matrix(0., (m*n,1)) ## f_ij corresponds to p[jm+i]
    I0 = set() #i in I0 <--> exists j s.t. P[i,j]==1
    J0 = set() #j in J0 <--> exists i s.t. P[i,j]==1
    for i in range(m):
        for j in range(n):
            p[j*m+i] = 1-P[i,j] # p[jm+i]==1: known path for [i,j] is absent
            if P[i,j]==1:
                I0.add(i); J0.add(j)

    mask = matrix(0., (m*n,1)) ## f_ij==0 corresponds to mask[jm+i]==1
    mask_indice = [] #list of (i,j) where mask_ij (mask[j*m+i]) is 1 (f_ij 0)
    for i in range(m):
        for j in range(n):
            if P[i,j]==0 and (i in I0 or j in J0):
                mask[j*m+i]=1
                mask_indice.append((i,j))

    useOldPD = True

    if len(mask_indice)>0:  #choice==1
        useOldPD = False

        #rhs = []
        #constraint from inedges
        I_w = []; J_w = [] #I/J indice constrained by edge weights w
        for i in range(m):
            prod = 1
            for j in range(n):
                prod = prod * mask[j*m+i]
            if prod==0:
                #rhs.append(a[i])
                I_w.append(i)
        #constraint from outedges
        for j in range(n):
            prod = 1
            for i in range(m):
                prod = prod * mask[j*m+i]
            if prod==0:
                #rhs.append(b[j])
                J_w.append(j)

        #prepare A and rhs
        A = matrix(0.0, (len(mask_indice), m*n))
        rhs = matrix(0.0, (len(mask_indice), 1))
        row_id = 0
        for i,j in mask_indice:
            A[row_id, j*m+i] = 1.0
            row_id += 1
        #pdb.set_trace()

        #prepare G and h
        #G = matrix(0.0, (m*n-len(mask_indice)+len(I_w)+len(J_w), m*n))
        #h = matrix(0.0, (m*n-len(mask_indice)+len(I_w)+len(J_w), 1))

        G = []
        h = []
        exclude_list = []

        row_id = 0

        for i in I_w: # sum_j f_ij >= w_i
            row = [0]*(m*n)
            row_indice = []
            for j in range(n):
                if  1: #(i,j) not in mask_indice:
                    row[j*m+i]=1.0
                    row_indice.append((i,j))
            G.append(row)
            h.append(float(a[i]))
            row_id += 1
            if len(row_indice)==1: #f_ij >= a[i], no need to constrain it >= 0 later
                exclude_list += row_indice
        #pdb.set_trace()

        for j in J_w: # sum_i f_ij >= w_j
            row = [0]*(m*n)
            row_indice = []
            for i in range(m):
                if 1: #(i,j) not in mask_indice:
                    row[j*m+i]=1.0
                    row_indice.append((i,j))
            G.append(row)
            h.append(float(b[j]))
            row_id += 1
            if len(row_indice)==1:
                exclude_list += row_indice
        #pdb.set_trace()
        
        for i in range(m):
            for j in range(n):
                '''if (i,j) in mask_indice or (i,j) in exclude_list:
                    continue ''' #f_ij=0, no need to further constrain here

                row = [0]*(m*n)
                row[j*m+i]=-1.0
                G.append(row)
                h.append(0.0)
                row_id += 1
        
        if (i,j) in mask_indice:
            row = [0]*(m*n)
            row[j*m+i]=1.0
            G.append(row)
            h.append(0.0)
            row_id += 1
        #pdb.set_trace()

        G = matrix(G).T
        h = matrix(h)
        #pdb.set_trace()       

        scale_factor = 1
        removal_factor = 0.4 #0.1;

        #if check_AbGh==1:
        #    pdb.set_trace()

    else:

        A = matrix(0.0, (m+n-1, m*n))
        rhs = matrix(0.0, (m+n-1, 1))

        row_id = 0
        for i in range(m):
            for j in range(n):
                A[row_id, j*m+i]=1.0
            rhs[row_id]=a[i]
            row_id += 1
        for j in range(n-1):
            for i in range(m):
                A[row_id, j*m+i]=1.0
            rhs[row_id]=b[j]
            row_id += 1

        G = spmatrix(-1.,range(m*n),range(m*n))
        h = matrix(0.,(m*n,1))

        scale_factor = max(max(rhs),1e-100) * 0.01;
        removal_factor = 0.4;

        #pdb.set_trace()

    ##########

    weight = LA.norm(a,1) ## (Not used)
    eps = 0.001
    tol = eps*weight  #test for significance.  Used for various purposes
    sparsity_factor = 0.4  #very aggressive curently - revert to 0.1 later
    #removal_factor = 0.4;
    scale = scale_factor #max(max(rhs),1e-100) * 0.01;
    #print(rhs)
    trials = int(round(min(2*m*n*max(m,n),100))) ## Number of randomized trials.
    curr_min = m*n +1  ## The lowest number of non known paths seen so far.
    curr_ans = [];     ## The best solution seen so far.
    curr_err = 0;
    curr_mult = 0;         ## multiplicity of current solution
    curr_on_unknown = 0;   ## amount of flow on non "known paths".
    for ctr in range(trials):  ## randomize the the coefficients for the known non=known path flow values.
        c=matrix(abs(numpy.random.normal(0,1,(m*n,1))))
        for i in range(m*n):
            if useOldPD == True:
                c[i] = c[i]*p[i]
            else:
                c[i] = -c[i]

        ## G and h are used to make sure that the flow values are non-negative.  
        #G = spmatrix(-1.,range(m*n),range(m*n))
        #h = matrix(0.,(m*n,1))
        if use_GLPK: 
            sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale,solver='glpk')
        else:
            if useOldPD == True:
                sol = solvers.lp(c =c ,G=G,h=h,A=A,b=rhs/scale)
            else:
                sol = solvers.lp(c =c ,G=G,h=h)
        temp_sol = array(sol['x'])*scale
        another_sol = copy.deepcopy(temp_sol)
        #print('useOldPD=%s temp_sol=\n%s'%(str(useOldPD), str(temp_sol)))
        if overwrite_norm: ## Do not use right now because a_true and b_true are wrong.
            #the true values of copy count are used to decide thresholding behavior
            a = a_true[:]
            b = b_true[:]
            
        ## This loop basically sets the flow values equal to 0 under certain conditions.              
        for i in range(m):
            for j in range(n):
                if another_sol[j*m+i]<sparsity_factor*min(a[i],b[j]):
                    another_sol[j*m+i]=0
                    
                if temp_sol[j*m+i]<removal_factor*min(a[i],b[j]) or temp_sol[j*m+i]<tol: # temporarily disabled
                    #if temp_sol[j*m+i] > 0.00001:
                    #    pdb.set_trace()
                    temp_sol[j*m+i]=0
                    another_sol[j*m+i]=0

                if another_sol[j*m+i]<0:
                    another_sol[j*m+i]=0
                    temp_sol[j*m+i]=0
        #print('useOldPD=%s temp_sol(after thresholding)=\n%s'%(str(useOldPD), str(temp_sol)))
        #pdb.set_trace()
        
        
        s=0 ## s equals how many non-zero flows we are sending down non "known paths" in the temp solution
        for i in range(m*n):
            if p[i] > 0:  #Only couont the paths that are not supported
                s=s+numpy.count_nonzero(another_sol[i])

        ## if the temporary solution is less than the current minimum solution, replace current minimum solution with
        ## the temporary solution.
        if s<curr_min:
            curr_min = s  ## current minimum value of non "known paths" in solution 
            curr_ans = temp_sol[:]  ## answer that attains it.  
            curr_mult = 0  ## This is how many times a solution with this many non-zero non-known paths is seen.
            curr_on_unknown = local_dot(array(p),temp_sol)  ## This says how much flow we are sending down non "known paths"

        else:
            if s==curr_min:
                if LA.norm(curr_ans-temp_sol) > tol:  ## Determines whethee to classify the solutions as different.
                    curr_mult = curr_mult +1
                if curr_ans==[] or ( abs(sum(sum(temp_sol))-sum(sum(curr_ans)))<tol and (local_dot(array(p),temp_sol) < curr_on_unknown) ) or sum(sum(temp_sol)) > sum(sum(curr_ans)):
                    ## These are a few conditions that make it the temporary solution:
                    ## 1: curr_sol is empty.  2: total flow difference is below a threshold AND less flow is going down unknown paths.  3: total flow is greater than curr_ans total flow.
                    curr_ans = temp_sol[:]
                    curr_on_unknown = local_dot(array(p),temp_sol)

    answer = matrix(0.,(m,n))
    if len(curr_ans) < m*n:
        pdb.set_trace()
    for i in range(m):
        for j in range(n):
            answer[i,j]=float(curr_ans[j*m+i])
    
    # Uniqueness consideration
    answer = array(answer)
    non_unique = 0
    if curr_mult > 1:
        non_unique =1
    
    # This makes the flow solution more sparse
    if sparsity != False:
        if m*n > sparsity:    
            tmp_dict = {}   
            for i in range(m):
                for j in range(n):
                    tmp_dict[(i,j)]=answer[i, j]
            sorted_tmp = sorted(tmp_dict.items(), key=operator.itemgetter(1))[::-1]
            sorted_tmp = sorted_tmp[:sparsity]
        
            new_ans = zeros((m, n))
            for ind in sorted_tmp:
                new_ans[ind[0][0], ind[0][1]] = ind[1]
            answer = new_ans 
    return [answer,non_unique]


#'''#Test Case
if __name__ == '__main__':

    import numpy;

    #choice = 0 # 0 old sf, 1 new sf

    '''f = [[1.5, 2],
         [3,   0]]

    P = [[0,   0],
         [1,   0]]'''

    '''f = [[0,   20,  0],
         [0,   30,  0],
         [15,  40,  21]]

    P = [[1,   0,  0],
         [1,   0,  0],
         [0,  0,  0]]'''

    '''
    #pdb.set_trace()

    check_AbGh = 0'''

    '''f = numpy.array(f)
    P = numpy.array(P)

    a = numpy.sum(f, axis=1)
    b = numpy.sum(f, axis=0)'''

    #exception case
    '''P = [[ 1.00e+00,  0.00e+00,  0.00e+00], 
         [ 0.00e+00,  1.00e+00,  0.00e+00],
         [ 0.00e+00,  0.00e+00,  0.00e+00]]
    P = numpy.array(P)
    a = [5.0, 179.99999972937738, 266.32720062622786]
    b = [19.970230104230318, 429.3599472409519, 1.9970230104230318]'''


    '''[a2,z] = path_decompose4(a,b,a,b,False,P,False,10)'''

    '''
    [m,n]=a2.shape
    newP = copy.deepcopy(a2)
    for i in range(m):
        for j in range(n):
        newP[i][j] = 1 if a2[i][j]>0 else 0
    '''
    #[b2,y] = path_decompose(a,b,a,b,False,newP,False,10)

    '''print(z)
    print('res=\n%s'%str(a2))'''
