import os
from datetime import datetime
import subprocess
import numpy as np
import sys

degree_folder = "../autogen/"+sys.argv[1].split('/')[-1]

def log(s):
    print(datetime.now(),"-",s,flush=True)

def execute(c):
        log("executing \""+c+"\"")
        s = subprocess.run(c,capture_output=True,shell=True)
        log("stdout result:\n"+s.stdout.decode('utf-8'))

def make(coulomb_radius,r0Be,gBe,ExBe,JC):
    name = "auto_"+coulomb_radius+"_"+r0Be+"_"+gBe+"_"+ExBe+"_"+JC
    exists = os.path.isdir(degree_folder+"/figures/"+name)
    if exists:
        log(name+" already exists, no reason to make again!")
    else:
        log("making "+name)
        execute("make -j8 autogen coulomb-radius="+coulomb_radius+" r0Be="+r0Be+" gBe="+gBe+" ExBe="+ExBe+" JC="+JC)
    return name

ExBe = []
gBe = []

for line in sys.stdin:
    line = line.rstrip()
    line = line.split()
    gBe.append(line[0])
    ExBe.append(line[1])
for Ex,g in zip(ExBe,gBe):
    a = float(Ex)
    b = float(g)
    #if a>2850 and a<2950 and b>34 and b<39:
    #make("12","1.42",g,Ex,"1-")
    make("8","1.42",g,Ex,"2-")
    #else:
        #log("skipping "+Ex+" "+g)

#if True:
    #for coulomb in ["8"]:
        #make(coulomb,"1.42","32.79","3037","2-")


