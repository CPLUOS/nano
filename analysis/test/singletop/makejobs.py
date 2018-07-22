#import sys, opt
import sys, os


# Usage : python makejobs.py [DEST DIR] j`python batch_nanoAOD.py n [# OF FILES PER JOB]`


strSubmitCardTemplate = """# Job description file for condor job samples_ZMM_PU140_pre4_fixed01_test01
executable = run_nanoAOD_to_ntuple.sh
universe   = vanilla
arguments  = $(Process) %(numPerJob)s %(dest)s
requirements = OpSysMajorVer == 6

log = %(jobname)s/condor.log

getenv     = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = %(jobname)s/job_$(Process).log
error = %(jobname)s/job_$(Process).err
transfer_input_files = job.tar.gz
#transfer_output_files =  singletop_nanoAOD/src/nano/analysis/test/singletop/res.root
#transfer_output_remaps = "res.root=%(jobname)s/res_$(Process).root"
queue %(jobnum)s
"""


if len(sys.argv) < 2 or sys.argv[ 1 ] == "--help":
  print "Usage : python makejobs.py [DEST DIR] j`python batch_nanoAOD.py n [# OF FILES PER JOB]`"
  sys.exit(0)

strJobNum = sys.argv[ 1 ]
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
strSubmitCard = strSubmitCardTemplate%{"dest": strJobNum, "jobname": strJobName, 
  "jobnum": strNumJob, "numPerJob": strNumRoot}
os.system("./full_archive.sh singletop_nanoAOD ; echo '@@ Archiving is complete'")
os.system("printf '%(card)s' | condor_submit"%{"card": strSubmitCard})


