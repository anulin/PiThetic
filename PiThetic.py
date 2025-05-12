#!/usr/bin/env python3
import subprocess
import sys
from itertools import combinations
from multiprocessing import Pool
import numpy as np
import scipy.optimize as op
from numpy import random as ra
from scipy.special import digamma, polygamma
from scipy.stats import poisson
import matplotlib.pyplot as plt
# from os import environ
# print(environ['OMP_NUM_THREADS'] )

def harmonic(n):
    return digamma(n+1)+np.euler_gamma
    return sum(1/(i) for i in range(1,n+1))

def htrzglik(a,b,p,ph):
    if b==3:
        p2=1-sum(p)
    else:
        p2=p[b]
    c,d={0,1,2,3}-{a,b}
    return p[a]*p2*np.prod([(1+2*pr)/6 for pr in ph[a]+ph[b]]+[(1-pr)/3 for pr in ph[d]+ph[c]])
def hmzglik(a,p,ph):

    b,c,d= {0, 1, 2, 3}-{a}
    if a==3:
        p=1-sum(p)
    else:
        p=p[a]
    return p**2*np.prod(ph[a])*np.prod([(1-pr)/3 for pr in ph[b]+ph[c]+ph[d]])

def loglikelihood(p, m,  ph):
    return -np.sum(np.log([p*pr+(1-p)*(1-pr)/3 for pr in ph[:m]]))-np.sum(np.log([p*(1-pr)/3+(1-p)*(pr) for pr in ph[m:]]))
def lhoodDiploidB(p, x, y, phs):
    c,d={b'a',b't',b'c',b'g'}-{x, y}
    if p>1 or p<0:
        return 0
    if type(phs)==list:
        return -np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) + 2 * p*(1-p)*np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                         (1-p) ** 2 * np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs])
    else:
        return -(p ** 2 * np.prod(phs[x]) * np.prod([(1 - pr) / 3 for pr in phs[y] + phs[c] + phs[d]]) + 2 * p*(1-p)*np.prod([(1+2*pr)/6 for pr in phs[x]+phs[y]]+[(1-pr)/3 for pr in phs[d]+phs[c]]) +
                 (1-p) ** 2 * np.prod([(1-pr) / 3 for pr in phs[x] + phs[c] + phs[d]]) * np.prod(phs[y]))
def mean(p,phs,x,y):

    c,d={b'a',b't',b'c',b'g'}-{x, y}

    p11=np.array([p**2* np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]])/
                 ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                   p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                   (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
    p12=np.array([ p * (1 - p) *2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]])/
                  (p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                   p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                   (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])


    return sum(2*p11+p12)/2/len(phs)-p
def meanfrqs(phs):
    vals=[]
    for x, y in combinations({b'a', b't', b'c', b'g'}, 2):
        if mean(1e-10,phs, x, y)>0 and mean(1-1e-10,phs, x, y)>0 :
            p=1
        elif mean(1e-10,phs, x, y)*mean(1-1e-10,phs, x, y)<0:
            p=op.root_scalar( mean, args=(phs, x, y), bracket=[1e-10,1-1e-10]).root
        else:
            p=0
        vals.append(p)
    return vals
def allelprb(ps,phs):
    prs = []
    for i in range(len(ps)):
        x, y=list(combinations({b'a', b't', b'c', b'g'}, 2))[i]
        c, d = {b'a', b't', b'c', b'g'} - {x, y}
        p=ps[i]
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs]))
        return prs
def hrm(phs,ps):
    prs = []
    vals = []
    for i in range(len(ps)):
        x, y = list(combinations({b'a', b't', b'c', b'g'}, 2))[i]
        c, d = {b'a', b't', b'c', b'g'} - {x, y}
        p = ps[i]
        p11 = np.array([p * 2 * np.prod(ph[x] / 2) * np.prod([(1 - pr) / 6 for pr in ph[y] + ph[c] + ph[d]]) /
                       (p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                        p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in
                                                                                                ph[d] + ph[c]]) +
                        (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph
                       in phs])
        p00 = np.array([(1 - p) * 2 * np.prod([(1 - pr) / 6 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y] / 2) /
                       (p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                        p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in
                                                                                                ph[d] + ph[c]]) +
                        (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph
                       in phs])
        vals.append(harmonic(sum(2 - p11 - p00)-1))
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod(
            [(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for
                           ph in phs]))
    return np.dot(prs,vals)/sum(prs)
def Dp(phs,ps, muts):
    prs=[]
    vals=[]

    for i in range(len(ps)):
        x, y = list(combinations({b'a', b't', b'c', b'g'}, 2))[i]
        c, d = {b'a', b't', b'c', b'g'} - {x, y}
        p = ps[i]

        # p11=np.array([p*2* np.prod(ph[x]/2) * np.prod([(1 - pr) / 6 for pr in ph[y] + ph[c] + ph[d]])/
        #              ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
        #                p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
        #                (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        # p00=np.array([(1-p)*2*np.prod([(1-pr) / 6 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]/2)/
        #              ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
        #                p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
        #                (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs]))
        n=2*len(phs)

        var=(muts-1)/muts*((11 * n ** 2 - 7 * n + 6) / 9 / n / (n - 1)) / (harmonic(n - 1) ** 2  -polygamma(1,n)+np.pi**2/6 )+(n+1)/3/(n-1)/muts-1/harmonic(n-1)**2
        vals.append(var)

    return np.dot(prs,vals)/sum(prs)
def theta(phs):
    prs=[]
    vals=[]
    sqrs=[]
    for x, y in combinations({b'a', b't', b'c', b'g'}, 2):
        c,d={b'a',b't',b'c',b'g'}-{x, y}
        if mean(1e-10,phs, x, y)>0 and mean(1-1e-10,phs, x, y)>0 :
            p=1
        elif mean(1e-10,phs, x, y)*mean(1-1e-10,phs, x, y)<0:
            p=op.root_scalar( mean, args=(phs, x, y), bracket=[1e-10,1-1e-10]).root
        else:
            p=0
        # p11=np.array(p*2* np.prod(ph[x]/2) * np.prod([(1 - pr) / 6 for pr in ph[y] + ph[c] + ph[d]])/
        #              ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
        #                p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
        #                (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs)
        # p00=np.array((1-p)*2*np.prod([(1-pr) / 6 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]/2)/
        #              ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
        #                p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
        #                (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs)
        # sz=sum(2-p11-p00)
        p11=np.array([p**2* np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]])/
                     ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                       (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs ])
        p00=np.array([(1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])/
                     ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                       (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs  ])
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs  ]))
        vals.append((1-np.prod(p11)-np.prod(p00))/harmonic(2*len(phs)-1))
        sqrs.append(vals[-1]/harmonic(2*len(phs)-1))
    avg=np.dot(prs,vals)/sum(prs)

    avgsq=np.dot(prs,sqrs)/sum(prs)
    return avg#/(avgsq-avg**2),1/(avgsq-avg**2)
def mutprob(phs,ps):
    prs=[]
    vals=[]
    for i in range(len(ps)):
        x, y = list(combinations({b'a', b't', b'c', b'g'}, 2))[i]
        c, d = {b'a', b't', b'c', b'g'} - {x, y}
        p = ps[i]

        p11=np.array([p**2* np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]])/
                     ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                       (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        p00=np.array([(1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])/
                     ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                       (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs]))

        vals.append(1-np.prod(p11)-np.prod(p00))

    return np.dot(prs,vals)/sum(prs)


def loglhoodDiploid(p,   ph):
    return -np.prod([ sum(hmzglik(i,p,ph)for i in range(4))+2*sum(htrzglik(a,b,p,ph) for a,b in combinations(range(4),2))])
def loglhoodPiDiploidB( phs):
    func=lambda p: sum(lhoodDiploidB(p, x, y, phs) + lhoodDiploidB(p, y, x, phs) for x, y in combinations({b'a', b't', b'c', b'g'}, 2))
    return op.minimize_scalar(func,bounds=[0,0.5])['x']
    # return loglhoodDiploidB(p, x, y, phs)+loglhoodDiploidB(p, y, x, phs)
    # return -np.prod([ p**2*np.prod(ph[0])*np.prod([(1-pr)/3 for pr in ph[1]])+2*p*(1-p)*htrzglik(0,1,p,ph)+
    #                      (1-p)**2*np.prod([(1-pr)/3 for pr in ph[0]])*np.prod( ph[1])])

def PAvg(phs,ps, pidata=None,**kwargs):
    prs=[]
    vals=[]
    if pidata==None:
        pidata = kwargs.get('c')
    for i in range(len(ps)):
        x, y = list(combinations({b'a', b't', b'c', b'g'}, 2))[i]
        c, d = {b'a', b't', b'c', b'g'} - {x, y}
        p = ps[i]

        p11=np.array([p**2* np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]])/
                     ( p**2*np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p*(1-p)*2 *np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]]) +
                       (1-p)**2*np.prod([(1-pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        p12=np.array( [2 *p*(1-p)*np.prod([(1+2*pr)/6 for pr in ph[x]+ph[y]]+[(1-pr)/3 for pr in ph[d]+ph[c]])/
                      (p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                       p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                       (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y])) for ph in phs])
        prs.append(np.prod([p ** 2 * np.prod(ph[x]) * np.prod([(1 - pr) / 3 for pr in ph[y] + ph[c] + ph[d]]) +
                           p * (1 - p) * 2 * np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y]] + [(1 - pr) / 3 for pr in ph[d] + ph[c]]) +
                           (1 - p) ** 2 * np.prod([(1 - pr) / 3 for pr in ph[x] + ph[c] + ph[d]]) * np.prod(ph[y]) for ph in phs]))
        #Vr=1/(4*p11+p12-(2*p11+p12)**2 )**2
        X=(2*p11+p12)#/Vr ; 2*p11/(p**2-p**4+p11-p11**2)
        p_2=(sum(X)**2-sum(X**2-2*p11))#/(sum(Vr)**2-sum(Vr**2-2/(p**2-p**4+p11-p11**2)))/4#/len(p11)/(2*len(p11)-1)/2#X^2-X==p^2*(n-1)*n
        # p_3=p_2*sum(X)-2*X**2*(sum(X)-X)-2*p11*X#..
        # p_4=p_2**2+2*(sum(X**2)**2-sum(X**4))-4*sum((X*(sum(X)-X))**2) -sum((2*p11)**2)-sum(X*(sum(X)-X)*8*p11)
        # Vr_pi=p_2-2*p_3+p_4-(sum(X)-p_2)**2

        vals.append((sum(X)/2/len(p11)-p_2/2/len(p11)/(2*len(p11)-1)))
    if pidata==None:
        return np.dot(prs,vals)/sum(prs)
    else:
        sp=sum(prs)
        pi=np.dot(prs,vals)/sp
        return pi/(sigma(phs,ps,prs)/sp-pi**2)/np.sum([1/(sigma(phs2,ps,prs)/sp-pi**2) for phs2 in pidata])

def sigma(phs,ps,prs):
    vals=[]
    #approximate pi square
    x,y,c,d={b'a', b't', b'c', b'g'}
    # coefs={ph:np.array([(1-cf)/cf for cf in ph[x]+ph[y] + ph[c] + ph[d]]) for ph in phs}
    #coefsp=1-3/(3+coefs)
    coefs = [np.array([(1 - cf) / cf for cf in ph[x] + ph[y] + ph[c] + ph[d]]) for ph in phs]
    prd = [np.prod(ph[x] + ph[y] + ph[c] + ph[d]) for ph in phs]
    bprd = [np.prod([(1 - pr) / 3 for pr in ph[x] + ph[y] + ph[c] + ph[d]]) for ph in phs]
    pprd = [np.prod([(1 + 2 * pr) / 6 for pr in ph[x] + ph[y] + ph[c] + ph[d]]) for ph in phs]

    for i in range(len(ps)):
        p = ps[i]

        X2=np.array([(p**2*prd[ph]-(1-p)**2*bprd[ph])**2/(prd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph])+sum(#no error ref
                    (p**2*prd[ph]*coefs[ph]/3-(1-p)**2*bprd[ph]/coefs[ph]*3)**2/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph]/coefs[ph]*3)+#1error ref->alt
                    2*((p**2*prd[ph]*coefs[ph]/3-(1-p)**2*bprd[ph])**2/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph])) +(1-p)**2 * bprd[ph])) +# 1error wrong
                    (p**2*bprd[ph]-(1-p)**2*prd[ph])**2/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph])+#no error alt
                    (p**2*bprd[ph]/coefs[ph]*3-(1-p)**2*prd[ph]*coefs[ph]/3)**2/(bprd[ph]/coefs[ph]*3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph]*coefs[ph]/3)+#1error alt->ref
                    2*((p**2*bprd[ph]-(1-p)**2*prd[ph]*coefs[ph]/3)**2/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph]))+ (1-p)**2 * prd[ph]*coefs[ph]/3))) +4*p-1 for ph in range(len(phs))])# 1error wrong

        Xp=np.array([(p**4*prd[ph]**2-(1-p)**2*p**2*bprd[ph]*prd[ph])/(prd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph])+sum(#no error ref
                    (p**4*prd[ph]*coefs[ph]**2/9-(1-p)**2*p**2*bprd[ph]*prd[ph])/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph]/coefs[ph]*3)+#1error ref->alt
                    2*((p**4*prd[ph]*coefs[ph]**2/9-(1-p)**2*p**2*bprd[ph]*prd[ph]*coefs[ph]/3)/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph]))+(1-p)**2 * bprd[ph])) +# 1error wrong
                    (p**4*bprd[ph]**2-(1-p)**2*p**2*prd[ph]*bprd[ph])/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph])+#no error alt
                    (p**4*bprd[ph]**2/coefs[ph]**2*9-(1-p)**2*p**2*prd[ph]*bprd[ph])/(bprd[ph]/coefs[ph]*3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph]*coefs[ph]/3)+#1error alt->ref
                    2*((p**4*bprd[ph]**2-(1-p)**2*p**2*prd[ph]*bprd[ph]*coefs[ph]/3)/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph])) +(1-p)**2 * prd[ph]*coefs[ph]/3))) +p**2 for ph in range(len(phs))])# 1error wrong
        p2=np.array([p**4*prd[ph]**2/(prd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph])+sum(#no error ref
                    p**4*prd[ph]*coefs[ph]**2/9/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * bprd[ph]/coefs[ph]*3)+#1error ref->alt
                    2*p**4*prd[ph]*coefs[ph]**2/9/(prd[ph]*coefs[ph]/3*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph]))+(1-p)**2 * bprd[ph]) +# 1error wrong
                    p**4*bprd[ph]**2/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph])+#no error alt
                    p**4*bprd[ph]**2/coefs[ph]**2*9/(bprd[ph]/coefs[ph]*3*p**2+ p * (1 - p) * 2 * pprd[ph]+ (1-p)**2 * prd[ph]*coefs[ph]/3)+#1error alt->ref
                    2*p**4*bprd[ph]**2/(bprd[ph]*p**2+ p * (1 - p) * 2 * pprd[ph]*2*(1-3/(3+coefs[ph]))+(1-p)**2 * prd[ph]*coefs[ph]/3)) +p**2 for ph in range(len(phs))])# 1error wrong
        N=2*len(phs)

        vals.append((sum(X2)*(2*N-1-4*p*(N-1)+4*p**2*(N-1)*(N-2)/(2*N-1))+4*sum(Xp)*(-1+2*p*(N-1)/(2*N-1))-8*p**3*N*(N-1)*(2*N-3))/N**2/(2*N-1)+
                    (N-1)/N*4*p**2+(sum(X2)**2-sum(X2**2))/N**2/(2*N-1)**2+(4*p**4*(N-1)*(4*(N-2)**2+1)+4*sum(p2)/N)/N/(2*N-1)**2)#had -pi**2 ??
    return np.dot(prs, vals)# / sum(prs)
def thetafindexp (nu,k,timesums):
    #np.mean(k*poisson.pmf(k,nu*timesums)-poisson.pmf(k+1,nu*timesums)*(k+1) )
    return np.mean(k*poisson.pmf(k,nu*timesums))/np.mean(poisson.pmf(k+1,nu*timesums)*(k+1))-1
def thetamedunb (nu,k,timesums):
    return np.mean(poisson.cdf(k,nu*timesums) )-0.5
# def mvarlglhood(arr,m1,m2,m3,ph):
#     p1,p2,p3=arr
#     return -np.sum(np.log([p1*pr+(1-p1)*(1-pr)/3 for pr in ph[:m1]]))-np.sum(np.log([p2*pr+(1-p2)*(1-pr)/3 for pr in ph[m1:m2]])) \
#            -np.sum(np.log([p3*pr+(1-p3)*(1-pr)/3 for pr in ph[m2:m3]]))-np.sum(np.log([(p2+p1+p3)*(1-pr)/3+(1-p1-p2-p3)*(pr) for pr in ph[m3:]]))
if __name__ == '__main__':
    if '--h' in sys.argv:
        print('''PiThetic.py --[flag]  [samtools mpileup input] 
        --h - help
        mode flags: 
            --freq - frequency (default)
            --pi - nucleotide diversity
            --pi [window] - calculate nucleotide diversity in window (non-MLE)
            --theta [window] - Watterson theta caclulated within windows of provided size
            --D [window] -calculates D'
            --accurate - improves accuracy for window statistics (only) D' and Pi but increases computation time
            --t - number of threads (default 1)
        for samtools mpileup input you may call samtools mpileup --h
        ''')
        exit()
    thrds=1
    tlsize=0
    windowD=False
    pi=False
    ff=False
    window=False
    strt=1#data start column
    cpi=0
    wpi=0
    positions=[]
    accurate=False
    if  '--t' in sys.argv :
        strt += 2
        thrds=int(sys.argv[sys.argv.index('--t')+1])
    if thrds>1:
        pool=Pool(thrds)
        tlsize=max(8,thrds)
    if  '--accurate' in sys.argv :
        strt += 1
        accurate=True
    if '--pi' in sys.argv :
        if sys.argv[sys.argv.index('--pi')+1].isdigit():
            strt+=1
            wpi=int(sys.argv[sys.argv.index('--pi')+1])

        strt+=1
        pi=True
    if '--theta' in sys.argv:
        window=int(sys.argv[sys.argv.index('--theta')+1])
        # step=int(sys.argv[sys.argv.index('--theta')+2])
        strt += 2
    if '--D' in sys.argv:
        windowD=int(sys.argv[sys.argv.index('--D')+1])
        # step=int(sys.argv[sys.argv.index('--theta')+2])
        strt += 2
    if not pi+window+windowD or '--freq' in sys.argv:
        if '--freq' in sys.argv:
            strt += 1
        ff=True
    arg=sys.argv[strt:]
    proc = subprocess.Popen(['samtools', 'mpileup',*arg],stdout=subprocess.PIPE)
    ccc=0
    muts=0
    chr=0
    thetas=0
    thetavrs=0
    pis=0
    times=[]
    timesums=[]
    N=10**6
    for i in range(2, sum(ii[-4:]=='.bam' for ii in sys.argv) + 1):
        times.append(ra.geometric(1 - (4*N-i*(i-1))/(4*N), size=200000))
    while True:
        reads=0
        stats = proc.stdout.readline()
        if not stats:
            break
        stats=stats.split(b'\t')
        phreds=[[1-10**(-0.1*(i-33)) for i in stats[j]] for j in range(5,len(stats),3)]
        nucleotides=stats[4::3]#4::3
        # print(stats[4::3])
        phs=[{b'a':[],b"g":[],b'c':[],b't':[]} for i in range(len(nucleotides))]
        ns=[{b'a':0,b"g":0,b'c':0} for i in range(len(nucleotides))]
        skip=0
        number=b''
        if window:
            ncld=b''
            #sampsz=0
            ccc += 1
        for i in range (len(nucleotides)):
            counter=0

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
                    # elif window and ncld!=nuc:
                    #     two=1
                    phs[i][nuc].append(phreds[i][counter])
                    counter+=1
                elif nuc.isalpha():
                    if window and ncld == b'':
                        ncld = nuc
                    # elif window and ncld != nuc:
                    #     two = 1
                    ns[i][nuc]+=1
                    phs[i][nuc].append(phreds[i][counter])
                    counter+=1
            reads+=counter
            # if window:
            #     if two:
            #         sampsz+=2
            #     else:
            #         sampsz+=2-2**-(len(phreds[i])-1)

        if wpi or windowD or tlsize:
            if chr!=stats[0]:
                pidata = []
                cpi=0
                ccc=1
                chr=stats[0]
            cpi += 1
            if reads>1:
                pidata.append(phs)
                if pi and not wpi:
                    positions.append((stats[0].decode(),stats[1].decode()))

        x,y=sorted( phs[0].keys(), key=lambda x:sum(len(phs[i][x]) for i in range(len(ns))) )[-2:]


        if window :
            if tlsize and ccc==window:
                thetas=sum(pool.map(theta,pidata))
                pidata=[]
            elif not tlsize and reads>1:
                tmp=theta(phs)
                thetas+= tmp#[0]
            #thetavrs+=tmp[1]
        # if int(stats[1]) >   b'23203': seem to be 0 reads
        #     print(nucleotides)
        if window and ccc == window:#"theta:",

            print(stats[0].decode(),str(int(stats[1])-ccc)+'-'+stats[1].decode(), thetas)#/thetavrs*window)
            ccc=0
            thetas=0
            thetavrs=0
        if ff:
            print(stats[0].decode(),stats[1].decode(), op.minimize_scalar(lhoodDiploidB, args=(y, x, phs), bounds=[0, 1])['x'], (y + b'(' + x + b')').decode())#'frequency:',
        if windowD:#"D':"
            # if not wpi and reads>1:
            #     pidata[-1].append(meanfrqs(phs))
            if windowD==cpi:
                if tlsize:
                    pidata=list(zip(pidata,pool.map(meanfrqs,pidata)))
                    muts=np.sum(pool.starmap(mutprob,pidata))
                    if muts != 0:
                        var=np.mean(pool.starmap(Dp,list(zip(*zip(*pidata),(muts,)*len(pidata)) )))
                else:
                    pidata=list(zip(pidata,map(meanfrqs,pidata)))
                    # hmc=np.mean([hrm(phs,ps) for phs,ps in pidata])
                    muts=np.sum([mutprob(phs,ps) for phs,ps in pidata])
                    #muts=sum(1/(1-muts))/sum(1/(muts*(1-muts)))*window
                    if muts!=0:
                        var=np.mean([Dp(phs,ps,muts)for phs,ps in pidata])
                if muts == 0:
                    print("no_mutations")
                    cpi = 0
                    pidata = []
                    continue
                if accurate:
                    Dpr=(sum(PAvg(phs,ps,c=pidata) for phs,ps in pidata)*window/muts-1/harmonic(2*len(phs)-1))/var
                else:
                    if tlsize:
                        Dpr = (np.sum(pool.starmap(PAvg,  pidata)) / muts - 1 / harmonic(2 * len(phs) - 1)) / var
                    else:
                        Dpr=(sum(PAvg(phs, p) for phs, p in pidata)/muts-1/harmonic(2*len(phs)-1))/var
                print(stats[0].decode(),str(int(stats[1])-windowD)+'-'+stats[1].decode(),Dpr)
                cpi=0
                pidata=[]
        if pi:
            if wpi  and cpi==wpi:
                pis+=1
                # pidata[-1].append(meanfrqs(phs))
                if reads>1:
                    print(stats[0].decode(),str(int(stats[1])-wpi)+'-'+stats[1].decode(),end=' ')
                    if not accurate:
                        if  tlsize:
                            print(np.mean(pool.starmap( PAvg, zip(pidata,pool.map(meanfrqs,pidata)) )) )#'nucleotide diversity:'
                        else:

                            print(np.mean([PAvg(phs, ps) for phs, ps in zip(pidata,map(meanfrqs,pidata))]))
                    elif tlsize:
                        print(np.sum(pool.starmap(PAvg , zip(pidata,pool.map(meanfrqs,pidata),(pidata,)*len(pidata)) )))
                    else:

                        print(sum(PAvg(phs,ps,c=pidata) for phs,ps in zip(pidata,map(meanfrqs,pidata))))#'nucleotide diversity:'
                cpi=0
                pidata=[]
            elif not wpi:#'nucleotide diversity:',
                if not tlsize:
                    if reads>1:
                        print(stats[0].decode(),stats[1].decode(),loglhoodPiDiploidB(phs))
                elif tlsize==cpi:
                    vls=pool.map(loglhoodPiDiploidB, pidata)
                    for vl in range(len(vls)): print(*positions[vl],vls[vl])
                    pidata,cpi,positions=[],0,[]



    if pi and not wpi and cpi!=0:
  # residual'nucleotide diversity:',
        vls = pool.map(loglhoodPiDiploidB,pidata)
        for vl in range(len(vls)): print(*positions[vl],vls[vl])
        pidata, cpi = [], 0

