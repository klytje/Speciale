import matplotlib.pyplot as plt
import sys
import math
from ROOT import TChain
from itertools import permutations
from dataframe import *

chain = TChain('a')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

r = dataframe(chain)
r = r.define("esum","eCM[0]+eCM[1]+eCM[2]")

xlim = (0,math.cos(math.pi/6))
ylim = (0,1.0)

plt.title(sys.argv[1])
r = r.filter("mul==3")\
    .filter("pT<35e3")\
    .filter("abs(deltaE)<200")\
    .define("E2","min(eCM[0],min(eCM[1],eCM[2]))")\
    .define("E1","max(min(eCM[0],eCM[1]), min(max(eCM[0],eCM[1]),eCM[2]))")\
    .define("E0","max(eCM[0],max(eCM[1],eCM[2]))")\
    .define("X","((E1-E2)/esum)*sqrt(3)")\
    .define("Y","3*E0/esum-1.0")\
    .filter("pow(X,2)+pow(Y,2)<1.0")\
    .filter("Y<0.93")\
    .filter("E2>250")

    #.filter("pT<50e3")\
    #.filter("abs(deltaE)<200")\
    #.define("E2","min(eCM[0],min(eCM[1],eCM[2]))")\
    #.define("E1","max(min(eCM[0],eCM[1]), min(max(eCM[0],eCM[1]),eCM[2]))")\
    #.define("E0","max(eCM[0],max(eCM[1],eCM[2]))")\
    #.define("X","((E0+2*E1)/esum-1)/sqrt(3)")\
    #.define("Y","E0/esum-1.0/3.0")\
    #.filter("pow(X,2)+pow(Y,2)<1.0/9.0")\
    #.filter("Y<0.31")
    #.filter("eCM[0]<6000 && eCM[1]<6000 && eCM[2]<6000")
#if "dt" in r.getColumnNames():
    #r = r.filter("dt<10e3 || dt>1e6")

ghist2d(r,"X","Y",100,(xlim,ylim),fname)
#if "wU" in r.getColumnNames():
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"f1"),weights="wU*f[0][0]")
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"f2"),weights="wU*f[0][1]")
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"im"),weights="wU*f[0][3]")
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"re"),weights="wU*f[0][2]")

#if "data" in fname:
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"count"))
#else:
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),export=(fname,"count"),weights="w0*w1*w2")
#
#plt.clf()
#if "wU" in r.getColumnNames():
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim),weights="wU*(0.22*f[0][0]+(1-0.22)*f[0][1]+2*sqrt(0.22*(1-0.22))*(f[0][2]*cos(2*3.1415*0.69)+f[0][3]*sin(2*3.1415*0.69)))*w0*w1*w2")
#else:
    #r.hist2d("X","Y",bins=100,range=(xlim,ylim))

plt.xlabel("$X$")
plt.ylabel("$Y$")
plt.gca().set_xlim(*xlim)
plt.gca().set_ylim(*ylim)
plt.gca().set_aspect("equal")

plt.savefig(fname+".pdf",dpi=600)

