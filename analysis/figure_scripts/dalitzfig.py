import matplotlib.pyplot as plt
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import *

chain = TChain('a')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

X = []
Y = []
for i,j,k in list(permutations([0,1,2])):
    X.append("((eCM["+str(j)+"]-eCM["+str(k)+"])/esum)*sqrt(3)")
    Y.append("3*eCM["+str(i)+"]/esum-1.0")
r = dataframe(chain)
r = r.define("esum","eCM[0]+eCM[1]+eCM[2]")

xlim = (-1.3,1.3)
ylim = (-1.3,1.3)

#plt.figure(figsize=(7*7/10,4.5*7/10))
width =  1.0/72.27*393
plt.figure(figsize=(width/2,width/1.61803398875*0.7))
##plt.title()
r = r.filter("mul==3")\
    .filter("pT<50e3")\
    .filter("abs(deltaE)<200")\
    .define("E2","min(eCM[0],min(eCM[1],eCM[2]))")\
    .define("E1","max(min(eCM[0],eCM[1]), min(max(eCM[0],eCM[1]),eCM[2]))")\
    .define("E0","max(eCM[0],max(eCM[1],eCM[2]))")\
    .define("X","((E1-E2)/esum)*sqrt(3)")\
    .define("Y","3*E0/esum-1.0")\
    .filter("pow(X,2)+pow(Y,2)<1.3")\
    .filter("E2>50")
    #.filter("Y<0.93")\

ghist2d(r,X,Y,200,(xlim,ylim),fname)
#if "wU" in r.getColumnNames():
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"f1"),weights=["wU*f[0][0]"]*6)
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"f2"),weights=["wU*f[0][1]"]*6)
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"im"),weights=["wU*f[0][3]"]*6)
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"re"),weights=["wU*f[0][2]"]*6)

#if "data" in fname:
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"count"))
    #print("data")
#else:
    #print("not data")
    #r\
    #.hist2d(X,Y,bins=200,range=(xlim,ylim),export=(fname,"count"),weights=["w0*w1*w2"]*6)

#plt.clf()
#if "wU" in r.getColumnNames():
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),weights=["wU*(0.22*f[0][0]+(1-0.22)*f[0][1]+2*sqrt(0.22*(1-0.22))*(f[0][2]*cos(2*3.1415*0.69)+f[0][3]*sin(2*3.1415*0.69)))*w0*w1*w2"]*6)
#elif "data" in fname:
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim))
#else:
    #r.hist2d(X,Y,bins=200,range=(xlim,ylim),weights=["w0*w1*w2"]*6)


#r.hist2d(X,Y,bins=250,range=(xlim,ylim),export=fname)
plt.xlabel("$X$")
plt.ylabel("$Y$")
plt.gca().zaxis.set_scale('log')
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)
plt.gca().set_aspect("equal")

plt.savefig(fname+".pdf",dpi=1200)

