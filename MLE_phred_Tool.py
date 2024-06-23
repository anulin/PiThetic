import subprocess
import sys
from scipy.stats import poisson
from numpy import sqrt
from itertools import combinations
import matplotlib.pyplot as plt
from numpy import random as ra
import scipy.optimize as op
import numpy as np
def htrzglik(a,b,p,ph):
    if b==3:
        p2=1-sum(p)
    else:
        p2=p[b]
    c,d={0,1,2,3}-{a,b}
    return p[a]*p2*np.prod([(1+2*pr)/6 for pr in ph[a]+ph[b]]+[(1-pr)/3 for pr in ph[d]+ph[c]])
def harmonic(n):
    return sum(1/(i) for i in range(1,n+1))
def hmzglik(a,p,ph):

    b,c,d= {0, 1, 2, 3}-{a}
    if a==3:
        p=1-sum(p)
    else:
        p=p[a]
    return p**2*np.prod(ph[a])*np.prod([(1-pr)/3 for pr in ph[b]+ph[c]+ph[d]])
def loglikelihood(p, m,  ph):
    return -np.sum(np.log([p*pr+(1-p)*(1-pr)/3 for pr in ph[:m]]))-np.sum(np.log([p*(1-pr)/3+(1-p)*(pr) for pr in ph[m:]]))
def loglhoodDiploidB(p, x, y, phs):
    c,d={b'a',b't',b'c',b'g'}-{x, y}
    if type(phs)==list:
        return -np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) + 2 * p*(1-p)*np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                        (1-p) ** 2 * np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs])
    else:
        return -(p ** 2 * np.prod(phs[x]) * np.prod([(1 - pr) / 3 for pr in phs[y] + phs[c] + phs[d]]) + 2 * p*(1-p)*np.prod([(1+2*pr)/6 for pr in phs[x]+phs[y]]+[(1-pr)/3 for pr in phs[d]+phs[c]]) +
                        (1-p) ** 2 * np.prod([(1-pr) / 3 for pr in phs[x] + phs[c] + phs[d]]) * np.prod(phs[y]))
def loglhoodDiploid(p,   ph):
    return -np.prod([ sum(hmzglik(i,p,ph)for i in range(4))+2*sum(htrzglik(a,b,p,ph) for a,b in combinations(range(4),2))])
def loglhoodPiDiploidB(p,x, y, phs):
    return sum(loglhoodDiploidB(p, x, y, phs)+loglhoodDiploidB(p, y, x, phs) for x,y in combinations({b'a',b't',b'c',b'g'},2))
    # return loglhoodDiploidB(p, x, y, phs)+loglhoodDiploidB(p, y, x, phs)
    # return -np.prod([ p**2*np.prod(ph[0])*np.prod([(1-pr)/3 for pr in ph[1]])+2*p*(1-p)*htrzglik(0,1,p,ph)+
    #                      (1-p)**2*np.prod([(1-pr)/3 for pr in ph[0]])*np.prod( ph[1])])
def thetafindexp (nu,k,timesums):
    #np.mean(k*poisson.pmf(k,nu*timesums)-poisson.pmf(k+1,nu*timesums)*(k+1) )
    return np.mean(k*poisson.pmf(k,nu*timesums))/np.mean(poisson.pmf(k+1,nu*timesums)*(k+1))-1
def thetamedunb (nu,k,timesums):
    return np.mean(poisson.cdf(k,nu*timesums) )-0.5
# def mvarlglhood(arr,m1,m2,m3,ph):
#     p1,p2,p3=arr
#     return -np.sum(np.log([p1*pr+(1-p1)*(1-pr)/3 for pr in ph[:m1]]))-np.sum(np.log([p2*pr+(1-p2)*(1-pr)/3 for pr in ph[m1:m2]])) \
#            -np.sum(np.log([p3*pr+(1-p3)*(1-pr)/3 for pr in ph[m2:m3]]))-np.sum(np.log([(p2+p1+p3)*(1-pr)/3+(1-p1-p2-p3)*(pr) for pr in ph[m3:]]))

if '--h' in sys.argv:
    print('''MLE_phred_Tool scoresfile --[flag]  [samtools mpileup input] 
    --h - help
    mode flags: 
        --freq - frequency (default)
        --pi - nucleotide diversity
        --theta [window] - Watterson theta caclulated within windows of provided size
    for samtools mpileup input you may call samtools mpileup --h
    ''')
    exit()
pi=False
ff=False
strt=1
if '--pi' in sys.argv :
    strt+=1
    pi=True
if '--theta' in sys.argv:
    window=int(sys.argv[sys.argv.index('--theta')+1])
    # step=int(sys.argv[sys.argv.index('--theta')+2])
    strt += 2
if not pi+window or '--freq' in sys.argv:
    if '--freq' in sys.argv:
        strt += 1
    ff=True
arg=sys.argv[strt:]
proc = subprocess.Popen(['samtools', 'mpileup',*arg],stdout=subprocess.PIPE)
ccc=0
muts=0
times=[]
timesums=[]
N=10**6
for i in range(2, sum(ii[-4:]=='.bam' for ii in sys.argv) + 1):
    times.append(ra.geometric(1 - (4*N-i*(i-1))/(4*N), size=200000))
while True:

    tns=0
    stats = proc.stdout.readline()
    if not stats:
        break
    stats=stats.split(b'\t')
    phreds=[[1-10**(-0.1*(i-33)) for i in stats[j]] for j in range(5,len(stats),3)]
    nucleotides=stats[4::3]#4::3

    # print(stats)
    phs=[{b'a':[],b"g":[],b'c':[],b't':[]} for i in range(len(nucleotides))]
    ns=[{b'a':0,b"g":0,b'c':0} for i in range(len(nucleotides))]
    skip=0
    number=b''
    if window:
        ncld=b''
        sampsz=0
        ccc += 1
    for i in range (len(nucleotides)):
        counter=0
        if window:
            two=0
        for nuc in np.frombuffer(nucleotides[i].lower(),dtype='S1'):#np.frombuffer so that nuc is b'' instead of int

            if skip!=0:
                skip-=1
                continue
            if number !=b'' and nuc.isalpha():
                skip=int(number)-1
                number=b''
            elif nuc.isdigit():
                if number==b'':
                    print('!!!!',b'\t'.join(stats))
                    exit()
                number+=nuc
            elif nuc==b'^':
                skip=1
    
            elif nuc== b'+' or nuc==b'-':
                number+=b'0'
            elif nuc== b't':
                # tns+=1
                if window and ncld == b'':
                    ncld = nuc
                elif window and ncld!=nuc:
                    two=1
                phs[i][nuc].append(phreds[i][counter])
                counter+=1
            elif nuc.isalpha():
                if window and ncld == b'':
                    ncld = nuc
                elif window and ncld != nuc:
                    two = 1
                ns[i][nuc]+=1
                phs[i][nuc].append(phreds[i][counter])
                counter+=1
        if window:
            if two:
                sampsz+=2
            else:
                sampsz+=2-2**-(len(phreds[i])-1)

    # phreds=[phreds[j] for j in phs[b"a"]]+[phreds[j] for j in phs[b'g']]+[phreds[j] for j in phs[b"c"]]+[phreds[j] for j in phs[b't']]
    # phs={key:[phreds[j] for j in phs[key]] for key in phs}

    # afrec=op.minimize_scalar(loglikelihood,args=(ns[b"a"],phreds),bounds=[0,1], method='Bounded')['x']    #
    # print ("A:", afrecs)
    # if counter==0:
    #     continue

    # allfrec=op.minimize(loglhoodDiploid,[0.25,0.25,0.25],args=list(phs[-1].values()),bounds=[[0,1]], constraints=op.LinearConstraint(np.diag([1,1,1]),0,1))
    # print(' '.join(str(i).upper() for i in zip(ns[-1].keys(), allfrec['x'])))
    x,y=sorted( phs[0].keys(), key=lambda x:sum(len(phs[i][x]) for i in range(len(ns))) )[-2:]


    if window and round(sampsz)>1:
        muts+=(op.minimize_scalar(loglhoodDiploidB, args=(y, x, phs), bounds=[0, 1])['x']<=1-1/round(sampsz))/harmonic(round(sampsz)-1)
    if ccc == window:
        print("theta:", muts)
        ccc=0
        muts=0
    if ff:
        print('frequency:',op.minimize_scalar(loglhoodDiploidB,args=(y,x,phs),bounds=[0,1])['x'],str(y+b'('+x+b')'))
    if pi:
        p=op.minimize_scalar(loglhoodPiDiploidB,args=(y,x,phs),bounds=[0,0.5])['x']
        print('nucleotide diversity:',p*(1-p)*2)
