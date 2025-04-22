import math
import cmath
# import numdifftools as nd
import numpy as np
import scipy.optimize as op
import statistics as st
import matplotlib.pyplot as plt
from scipy.stats import binom
from scipy.misc import derivative as deriv
from decimal import *
from scipy.integrate import tplquad, nquad
from time import time
def binplus(n,k):#2n-polynomial order, k - order of x
    n-=1
    k-=1
    return 2**(k-2*n)*sum(math.comb(i,k)*math.comb(2*n,2*i) for i in range(k,n+1))

def lnlikepi(pi,n,r):
    return np.log((1+np.sqrt(1-2*pi))**(n-r)*(1-np.sqrt(1-2*pi))**r+(1-np.sqrt(1-2*pi))**(n-r)*(1+np.sqrt(1-2*pi))**r)
def lnp(pi1,a,r1,r2):

    return lnlikepi(pi1,20,r1)+lnlikepi(a+pi1,20,r2)
def lnneg(pi1,a,r1,r2):

    return lnlikepi(pi1,20,r1)+lnlikepi(a-pi1,20,r2)
def like(x,n,r,nozero=True):
    return -(binom.pmf(n - r, n, x) + binom.pmf(r, n, x))
xi=0.35
print(op.minimize_scalar(lambda pi: -(lnlikepi(pi*xi,20,3)+lnlikepi(pi*(1-xi),20,7)), method='Bounded',bounds=[0,1]))
print(op.minimize_scalar(lambda pi:-lnlikepi(pi,20,7), method='Bounded',bounds=[0,0.5])['x'])

plt.plot(np.linspace(0,1,1000),[lnlikepi(pi*xi,20,3)+lnlikepi(pi*(1-xi),20,7) for pi in np.linspace(0,1,1000)])
plt.show()

vals=[]
for xi in np.linspace(0,1,100):
    vals.append(op.minimize_scalar(lambda pi: -(lnlikepi(pi*xi,20,3)+lnlikepi(pi*(1-xi),20,7)), method='Bounded',bounds=[0,1])['x'])
plt.plot(np.linspace(0,1,100)*vals,(1-n.linspace(0,1,100)*vals))
plt.show()
print(vals)
exit()
#pi2=a-pi1
print(lnlikepi(0,20,3))

alp=0.01
f=[[0]]
for i in np.linspace(0,1,1000):
    a=(deriv(lnlikepi,alp,args=(20,3),dx=0.0001)-deriv(lnlikepi,alp,args=(20,7),dx=0.0001))
    a1=deriv(lnlikepi,alp,args=(20,3),dx=0.0001,n=2)+deriv(lnlikepi,alp,args=(20,7),dx=0.0001,n=2)#deriv(lnneg,alp,args=(alp,3,7),dx=0.0001,n=2)
    b=deriv(lnlikepi,alp,args=(20,3),dx=0.0001)-deriv(lnlikepi,alp,args=(20,7),dx=0.0001,n=2)

    f[-1].append( (a-a1)/a*f[-1][-1]-b/a)
def b_a(pi1, pi):
    return (deriv(lnlikepi,pi1,args=(20,3),dx=0.0001,n=2)-deriv(lnlikepi,pi-pi1,args=(20,7),dx=0.0001,n=2))/(deriv(lnlikepi,pi1,args=(20,3),dx=0.0001)-deriv(lnlikepi,pi-pi1,args=(20,7),dx=0.0001))#=b/a

for pi in np.linspace(0.6,0.601,7):
    print([b_a(pi1,pi) for pi1 in np.linspace(0,pi,1000)])
    plt.plot(np.linspace(0,pi,1000),[b_a(pi1,pi) for pi1 in np.linspace(0,pi,1000)])
plt.show()
print(deriv(lambda x:lnlikepi(x,20,3)+lnlikepi(a-x,20,7),a,args=(20,3), dx=0.0001,n=2))
print(np.log(np.array([1,2,3])))
exit()
MLV=[]
MlM=[]
SlM=[]
MLLM=[]
SLV=[]
ML=[]
SL=[]
n=20

def p_pi(pi):
    if pi>0.5:

        print(pi,)

    if pi<0:
        return 1
    return (1+cmath.sqrt(1-2*pi))/2

if n%2==1:
    coeffs=np.empty(((n+1)//2,(n+1)//2))
else:
    coeffs=np.empty((n//2,n//2))



coeffseven=np.array([1])*2
coeffsodd=np.array([2/3])*3
if n%2==0:

    coeffs[1]=np.round(np.append([0]*(len(coeffs)-1),coeffseven),2)
    i=1
    print(np.round(coeffs[i] * ([0] * (10 - i) + [1 / math.comb(2*i, j+1) for j in range(i)]), 2), binplus(i, 1), i)
    for i in range(2,n//2):

        coeffs[i]=np.sum([2*binplus(i,j+1)/binplus(i,1)*coeffs[j+1]*(-1)**(j+i) for j in range(i-1)], axis=0)#{c(n,k)}/c(n,1)(==n)*2**(n//2-1)

        coeffs[i,-i]=1/binplus(i,1)
        print(np.round(coeffs[i]*([0]*(10-i)+[1/math.comb(2*i,j+1) for j in range(i)]), 2), binplus(i, 1), i)
print(np.round(coeffs,2))
exit()

print(binom.pmf(81, 81, 0)+binom.pmf(0,81 , 0),like(0,81,0))

print(deriv(lambda x,a :a*x**3,1,0.0001,2,order=7, args=[11]))
#print(deriv(lambda x,n,r: like(x,n,r), 0.5,0.0001,2,args=(81,30)))
plt.plot(np.linspace(0,0.5,31),[like(x,20,0)+like(x,20,5)  for x in np.linspace(0,0.5,31)])#sum(like(x,20,n) for n in range(2,11))
plt.plot(np.linspace(0,0.5,31),[like(x,20,5) for x in np.linspace(0,0.5,31)])
plt.plot(np.linspace(0,0.5,31),[like(x,20,6) for x in np.linspace(0,0.5,31)])
plt.show()
p=op.minimize_scalar(like,0.38,args=(4,3, False),bounds=[0,0.5], method='Bounded')['x']#
print(p,2*(p-p**2))

plt.plot(np.linspace(0.4,0.5,55),list(zip([derivlogLpi(i,81,40,2)for i in np.linspace(0.4,0.5,55)],[derivlogLpi(i,81,37,2)for i in np.linspace(0.4,0.5,55)],[derivlogLpi(i,81,43,2)for i in np.linspace(0.4,0.5,55)])))
plt.show()
mlePis=[]
comps=[[]]
comps2=[[]]
vals=[[]]
n2=101
n1=81

rrs=[4,5,6,5,4,6]
depth=len(rrs)
def tst(ns,*args):
    print(sum(args)+ns)
tst(2,3,2,3,2,7)
tm=time()
# ttt=binom.rvs(40, 0.44, size=int(5e4))/40
# print(np.mean(2*(ttt-ttt**2)*40/39))
# print(np.mean([sum(x!=y for x,y in np.array_split(binom.rvs(1, 0.44, size=40),20)) for i in range(int(5e4))])/20)
# print(nquad(lambda *args:np.prod([likepi(args[i],args[-2],args[:-1][i]) for i in range(len(args[:-2]))]), [lambda *args: (min(0.3 -sum(args[:-2]),0.5), max(0.3-sum(args[:-2])-0.5 ,0))]*depth, args=(20,rrs))[0]/
#       nquad(lambda *args:1 , [lambda *args: (min(0.3 -sum(args[:-2]),0.5), max(0.3-sum(args[:-2])-0.5 ,0))]*depth, args=(20,rrs))[0])
# print(time()-tm)
# exit()

# frqs=[1-p_pi(2 * (i / n1 - i ** 2 / n1**2) *n1/(n1-1)).real for i in range(21)][:-3]
frqs=[ round(op.minimize_scalar(like, method='Bounded',args=(n1,i, False),bounds=[0,1])['x'],4) for i in range (n1//2+1)][:-3]
for n in [x**2 for x in range(9,10)]:#(2,10)
    precalc={}
    precalc2={}
    for prob in frqs:#(0.1,0.5,5) (0.3,0.3,1)np.linspace(0.01,0.25,12)np.linspace(0.0,0.5,11)
        ML=[]
        MLL=[]
        SL=[]
        rs=binom.rvs(n, prob, size=500000)
        ks=binom.rvs(n1, prob, size=500000)
        ks2=binom.rvs(n2, prob, size=500000)
        for i in range(50000):
            L=[]
            r=rs[i]#0.3
            k=ks[i]
            k2=ks2[i]
            pi_n= 2 * (k2 / n2 - (k2 / n2) ** 2) #*n2/(n2-1)
            if k in precalc2:
                pi100=precalc2[k]
            else:
                p=(op.minimize_scalar(like,method='Bounded',args=(n1,k, False),bounds=[0,1])['x'])
                # pi100= 2 * (k / n1 - k ** 2 / n1**2) *n1/(n1-1)#*100/99
                # if pi10>0.5:
                #     pi10=0.5

                pi100=2*(p-p**2)
                precalc2[k]= pi100

            comps[-1].append(pi100 < pi_n)
            comps2[-1].append(pi100 > pi_n)
            vals[-1].append(pi100)
            if r in precalc:
                p,b,bb=precalc[r]
            else:
                p=(op.fmin_tnc(like,r/n,args=(n,r, False),approx_grad=True,bounds=[(0,1)],disp=0)[0][0])
                pi=(2*(p-p**2))
                bb = bias[int((r+2) / 82 * 12)-1]
                #print(pi)
                likes=[-like(p,n,j)/2 for j in range(n+1)]

                # I=sum(-derivlogLpi(pi,n,j,2)*likes[j] for j in range(len(likes))).real
                #print(max(-derivlogLpi(pi,n,j,2)*likes[j] for j in range(len(likes))))
                # J=sum(derivlogLpi(pi,n,j,1)*derivlogLpi(pi,n,j,2)*likes[j] for j in range(len(likes))).real
                #K=sum(derivlogLpi(pi,n,j,3)*likes[j]  for j in range(len(likes))).real
                #     I=-derivlogL(p,n,r,2)
                #     J=derivlogL(p,n,r,1)
                #     K=derivlogL(p,n,r,3)
                #     print(I,J,K)
                #     if math.isnan(I):
                #         print(pi,p,n,r,likes[0])

                #b=0.5/(I**2)*(K+2*J)
                b=0
                if(bb>0.1):
                    print(bb,r)
                precalc[r]=(p,b,bb)

            #print(b,r,p)
            ML.append(1-(p)**2-(1-p)**2-b)
            MLL.append(1-(r/n)**2-(1-r/n)**2)
            if b>0.5:
                print(sum(likes))
                #   print(b,p,r,I,J,K)
                exit()
                #print([-derivlogL(p,n,j,2)*likes[j] for j in range(len(likes))])
            if 2*(r/n-n/(n-1)*((r**2-r)/n**2))>0.5:
                SL.append(0.5)
            else:
                SL.append(2*(r/n-n/(n-1)*((r**2-r)/n**2)))
        if n==81:
            MlM.append(st.mean(i-2*prob*(1-prob)for i in ML))
            MLLM.append(st.mean(i-2*prob*(1-prob) for i in MLL))
            SlM.append(st.mean(i for i in SL)-2*prob*(1-prob))
        print(st.variance(ML), st.variance(SL), st.variance(MLL))
        print(p,r/n,L)
        print(MlM[-1],SlM[-1],2*(r/n-(r/n)**2),MLLM[-1])
        print(2*(prob-prob**2))
        #exit()####!
        MLV.append([st.variance(ML),st.variance(SL),st.variance(MLL)])
        print('compare pi100<pin, >:',sum(comps[-1]), sum(comps2[-1]))
        comps.append([])
        comps2.append([])
        vals.append([])
    if n==81:
        print(MlM)

        # plt.plot( np.linspace(0.01,0.25,12),list(zip(MLLM,SlM)))
        # print(list(zip(MlM,MLV)))
        # plt.fill_between(np.linspace(0.01,0.25,12), [a_i - b_i[2] for a_i, b_i in zip(MLLM,MLV)], [a_i + b_i[2] for a_i, b_i in zip(MLLM,MLV)], alpha=0.2, label='error band')
        # plt.fill_between(np.linspace(0.01, 0.25, 12), [a_i - b_i[1] for a_i, b_i in zip(SlM, MLV)],
        #                  [a_i + b_i[1] for a_i, b_i in zip(SlM, MLV)], alpha=0.2, label='error band')
        # # plt.fill_between(np.linspace(0.1, 0.5, 12), [a_i - b_i[2] for a_i, b_i in zip(MLLM, MLV)],
        # #                  [a_i + b_i[2] for a_i, b_i in zip(MLLM, MLV)], alpha=0.2, label='error band')
        # plt.legend(['Classic', 'Corrected for df'])
        # plt.show()
        # plt.plot( np.linspace(0.01,0.25,12),list(zip(MlM,MLLM)))
        # print(list(zip(MlM,MLV)))
        # plt.fill_between(np.linspace(0.01,0.25,12), [a_i - b_i[0] for a_i, b_i in zip(MlM,MLV)], [a_i + b_i[0] for a_i, b_i in zip(MlM,MLV)], alpha=0.2, label='error band')
        # plt.fill_between(np.linspace(0.01, 0.25, 12), [a_i - b_i[2] for a_i, b_i in zip(MLLM, MLV)],
        #                  [a_i + b_i[2] for a_i, b_i in zip(MLLM, MLV)], alpha=0.2, label='error band')
        # # plt.fill_between(np.linspace(0.1, 0.5, 12), [a_i - b_i[2] for a_i, b_i in zip(MLLM, MLV)],
        # #                  [a_i + b_i[2] for a_i, b_i in zip(MLLM, MLV)], alpha=0.2, label='error band')
        # plt.legend(['Corrected MLE', 'Classic'])
        # plt.show()
#print(vals)
print(precalc2)
# print(vals[0])
bx=plt.boxplot(vals[:-1],labels=[round(i,2) for i in frqs],notch =True, medianprops={'linewidth':4, 'color':'r'})
sc=plt.plot(range(1,len(frqs)+1),[2*(i-i**2) for i in frqs],linewidth=1,marker='o',color='blue', label='true value')#380
plt.xlabel('Allele frequency')
plt.ylabel('π')
plt.legend([bx['medians'][0],sc[0]], ['median','true value'],loc='lower right')
plt.show()
pis=[2*(i-i**2)for i in np.linspace(0.01,0.25,12)]
diffs=[[abs(vals[i][j]-pis[i]) for j in range(len(vals[0]))]for i in range(12)]
plt.boxplot(diffs,labels=[round(i,2) for i in np.linspace(0.01,0.25,12)],notch =True)
#plt.scatter(range(1,13),[2*(i-i**2) for i in np.linspace(0.01,0.25,12)],linewidth=1,marker='_',color='green',s=380)#
plt.xlabel('p')
plt.ylabel('|П^-П|')
plt.show()
plt.plot(np.linspace(0.01,0.25,12),[[i[0],i[2]] for i in MLV])#(2,10)
plt.show()

