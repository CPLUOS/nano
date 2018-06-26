#import sys, opt
import sys, os


# Usage : python makejobs.py nanoAOD j`python batch_nanoAOD.py n [# OF FILES PER JOB]`


if len(sys.argv) < 2 or sys.argv[ 1 ] == "--help":
  print "Usage : python makejobs.py nanoAOD j`python batch_nanoAOD.py n [# OF FILES PER JOB]`"
  sys.exit(0)

strJobNum = sys.argv[ 1 ]

strJobName = "res_" + strJobNum
strJobName = "nanoAOD_to_ntuple"

#strNumJob = "500" if len(sys.argv) < 3 else sys.argv[ 2 ]
strNumJob  = "250"
strNumRoot = "2"

if len(sys.argv) > 2: 
  for i in range(len(sys.argv)): 
    if i == 0 or i == 1: continue
    #if sys.argv[ i ][ 0 ] == "j": strNumJob = sys.argv[ i ][ 1: ]
    #if sys.argv[ i ][ 0 ] == "n": strEvtNum = sys.argv[ i ][ 1: ]
    if sys.argv[ i ][ 0 ] == "j": 
      strNumRoot = sys.argv[ i ][ 1: ].split(",")[ 0 ]
      strNumJob  = sys.argv[ i ][ 1: ].split(",")[ 1 ]

if not os.path.isdir(strJobName): 
  os.mkdir(strJobName)
  open(strJobName + "/.create-batch", "w").write("")

#open("run.sh", "w").write(open("run_template.sh", "r").read()%{"runnumber": strJobNum, "samplesize": strEvtNum})
open("submit.jds", "w").write(open("submit_nanoAOD_to_ntuple_template.jds", "r").read()%{"jobname": strJobName, "jobnum": strNumJob, "numPerJob": strNumRoot})

os.system("./full_archive.sh singletop_nanoAOD ; echo '@@ Archiving is complete' ; condor_submit submit.jds"%{"runnumber": strJobNum, "strJobName": strJobName})


