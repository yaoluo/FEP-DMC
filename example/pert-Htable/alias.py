import numpy as np
def construct_alias(weight):
    N = len(weight)
    Q = weight/np.sum(weight) * N
    Q_c = Q+0
    print(N)
    index = np.arange(N)
    small = []
    large = []
    balanced = []
    alias = -np.ones(N).astype(int)
    
    for i,x in enumerate(Q):
        if x>1:
            large.append(i)
        elif x<1:
            small.append(i)
        else :
            alias[i] = i #alias of itself 
    
    while len(small)!=0 and len(large)!=0:
        sn = small.pop()
        ln = large.pop()
        
        alias[sn] = ln 
    
        Q[ln] = Q[ln] + Q[sn] - 1
        if Q[ln]>1:
            large.append(ln)
        elif Q[ln]<1:
            small.append(ln)
        else:
            alias[ln] = ln 
            #alias of itself 
    for sn in small:
        alias[sn] = sn 
    for ln in large:
        alias[ln] = ln 
    
    #check alias table 
    Q_inv = np.zeros(N)
    for i in range(N):
        Q_inv[i] = Q_inv[i] + Q[i]
        Q_inv[alias[i]] = Q_inv[alias[i]] + 1-Q[i]
    print('error = ',np.max(np.abs(np.asarray(Q_inv)-np.asarray(Q_c))))
    return alias,Q


import sys 
args=str(sys.argv)

dirs=sys.argv[1]
print(sys.argv[1])

x = np.loadtxt(dirs+'/qgrid-dmc.dat')
x[:,3] = x[:,3] / np.sum(x[:,3]) * len(x[:,3])
alias,Q = construct_alias(x[:,3])

file=open(dirs+'/alias.dat','w')
for i in range(len(Q)):
    #print(f'{alias[i]}       {Q[i]:20.10E}')
    file.write(f'{alias[i]+1}       {Q[i]:20.10E}   {x[i,3]:20.10E}\n')

file.close()
