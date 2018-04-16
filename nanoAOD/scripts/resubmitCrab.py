#!/usr/bin/env python
# --siteblacklist=T2_ES_IFCA,T2_US_Florida
import sys,os,time

moreopts=''
for o in sys.argv:
    if o == sys.argv[0]:
        continue
    moreopts+=' '+o

print moreopts

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
        runCommand("crab resubmit -d "+d+moreopts, n)
    os.system("hostname")
    time.sleep(5000)

