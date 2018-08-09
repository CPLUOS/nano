import os, sys, json


strPathDraw = "%s/src/nano/analysis/test/singletop/draw"%os.environ[ "CMSSW_BASE" ]
dicSetDef = json.load(open(os.path.join(strPathDraw, "listSet.json")))

strArg = "".join(s + " " for s in sys.argv[ 2: ])
strVars = " ".join("\"%s\""%s for s in dicSetDef[ "Vars" ].keys())

strCode = """
for vars in %(vars)s ; do 
  python plotDraw.py %(src)s -p \"$vars\" %(args)s & 
  while true ; do 
    if [ `ps -ef | grep $UID | grep python | grep plotDraw | wc -l` -lt 16 ] ; then break ; fi
    sleep 1
  done
done
"""%{"src": sys.argv[ 1 ], "vars": strVars, "args": strArg}

os.system("echo '" + strCode + "' | bash")

#for strVar in dicSetDef[ "Vars" ].keys(): 
#  os.system("python plotDraw.py %s -p \"%s\" %s"%(sys.argv[ 1 ], strVar, strArg))


