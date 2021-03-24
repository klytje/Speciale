import matplotlib.pyplot as plt
import sys
from ROOT import TChain
from dataframe import *

chain = TChain('tree')
for f in input().split(): chain.AddFile(f)
fname = sys.argv[1]+"/"+sys.argv[0].split('/')[2].split('.')[0]

r = dataframe(chain)
r = r.define("esum","E_cm[0]+E_cm[1]+E_cm[2]")
        #.define("X","((E_cm[0]+2*E_cm[1])/esum-1)/sqrt(3)")\
        #.define("Y","E_cm[0]/esum-1.0/3.0")\
        #.define("phi","fmod(atan2(Y,X)+TMath::Pi(),2*TMath::Pi()/3)")\
        #.filter("pow(X,2)+pow(Y,2)<1.0/9.0")\
lim = (0,3.1416/3)
r = r.filter("p_tot<35e3")\
    .filter("abs(deltaE)<200")\
    .define("E2","min(E_cm[0],min(E_cm[1],E_cm[2]))")\
    .define("E1","max(min(E_cm[0],E_cm[1]), min(max(E_cm[0],E_cm[1]),E_cm[2]))")\
    .define("E0","max(E_cm[0],max(E_cm[1],E_cm[2]))")\
    .define("X","((E1-E2)/esum)*sqrt(3)")\
    .define("Y","3*E0/esum-1.0")\
    .filter("pow(X,2)+pow(Y,2)<1.0")\
    .filter("Y<0.93")\
    .filter("E2>250")

ghist(r,"atan2(X,Y)",300,lim,fname)
#if "wU" in r.getColumnNames():
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"f1"),weights="wU*f[0][0]")
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"f2"),weights="wU*f[0][1]")
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"im"),weights="wU*f[0][3]")
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"re"),weights="wU*f[0][2]")

#if "data" in fname:
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"count"))
#else:
    #r.hist("atan2(X,Y)",bins=300,range=(lim),export=(fname,"count"),weights="w0*w1*w2")
#plt.clf()
#if "wU" in r.getColumnNames():
    #r.hist("atan2(X,Y)",bins=300,range=(lim),weights="wU*(0.22*f[0][0]+(1-0.22)*f[0][1]+2*sqrt(0.22*(1-0.22))*(f[0][2]*cos(2*3.1415*0.69)+f[0][3]*sin(2*3.1415*0.69)))*w0*w1*w2")
#else:
    #r.hist("atan2(X,Y)",bins=300,range=(lim))

plt.title(sys.argv[1])
plt.xlabel("$\\phi$")
plt.ylabel("Count")

plt.savefig(fname+".pdf",dpi=600)
