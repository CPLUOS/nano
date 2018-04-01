#!/usr/bin/env python

import os,time

def runCommand(commandstring, loopn):
    print "loop", loopn, ",", commandstring
    os.system(commandstring)
    
multicrabdir = []
for fn in os.listdir('.'):
    if os.path.isdir(fn):
        if fn.startswith('crab_'):
            multicrabdir.append(fn)
            print "found crab3 job dir:", fn

for n in range(20):
    for d in multicrabdir:
        runCommand("crab resubmit -d "+d, n)
    os.system("hostname")
    time.sleep(5000)

