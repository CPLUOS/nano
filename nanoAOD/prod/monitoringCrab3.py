#!/usr/bin/env python
from CRABAPI.RawCommand import crabCommand
import glob,sys,optparse,copy,time,httplib

class CrabMonitor :
  def __init__(self, argument) :
    self.jobs=glob.glob("%s*"%argument)
    self.remain_jobs = copy.deepcopy(self.jobs)
    self.completed_jobs = []
    self.submitfailed_jobs=[]
    self.jobs_status ={}
    self.statuses ={}
    if ( len(self.jobs) ==0 ) :
      print "no jobs"
      sys.exit(-1)

  def cleaningJob(self) :
    completed_jobs = [job for job in self.remain_jobs if (self.statuses.get(job,"No info") == "COMPLETED")]
    submitfailed_jobs = [job for job in self.remain_jobs if (self.statuses.get(job,"No info") == "SUBMITFAILED")]
    remain_jobs = [job for job in self.remain_jobs if (job not in completed_jobs and job not in submitfailed_jobs)]
    self.remain_jobs = remain_jobs
    self.completed_jobs = completed_jobs
    self.submitfailed_Jobs = submitfailed_jobs

  def status(self ) :
    for idx, job in enumerate(self.remain_jobs) :
      print "Jobs (%d/%d)"%(idx+1,len(self.remain_jobs))
      while(True) :
        try : 
          res = crabCommand('status', job)
        except httplib.HTTPException, e:
          print "HTTPException",e," try again!"
          continue
        self.statuses[job]=res['status']
        self.jobs_status[job]=res['jobsPerStatus']
        print self.statuses[job], self.jobs_status[job]
        if ( self.statuses.get(job,"No info") != "No info" and self.jobs_status.get(job,"No info")!="No info" ) :
          break
    self.cleaningJob()
   
  def printAllStatus(self):
    print "All jobs's status : "
    for job in self.jobs :
      status = self.statuses.get(job,"No info")
      print job,status 
  def printAllJobsStatus(self) :
    print "All jobs's task status : "
    for job in self.jobs :
      status = self.jobs_status.get(job,"No info")
      print job,status
  def printRemainJobsStatus(self) :
    print "Remain jobs's task status : "
    for job in self.remain_jobs :
      status1 = self.statuses.get(job,"No info")
      status2 = self.jobs_status.get(job,"No info")
      print job,status1, status2
  def printCompletedJobs(self) :
    print "Complited jobs : "
    for job in self.completed_jobs :
      print job
  def totalNumOfStatus(self) :
    totalStatus={}
    for job in self.remain_jobs :
      status = self.jobs_status.get(job,"No info")
      if ( status == "No info" ) :
        continue
      for st in status.keys() :
        if ( totalStatus.get(st,"NoInfo") != "NoInfo") :
          totalStatus[st] += status[st]
        else :
          totalStatus[st] = status[st]
    return totalStatus
  def collectFailed(self) :
    status = copy.deepcopy(self.jobs_status)
    sort_failed = sorted( status, key=status.__getitem__)
    return sort_failed
  def Monitoring(self,maxJob) :
    while( True ) :
      if ( len( self.remain_jobs) == 0 ) :
        break
      print self.submitfailed_jobs
      self.status()
      t_status = self.totalNumOfStatus()
      print t_status
      num_failed = t_status.get('failed',0)
      num_transferring = t_status.get('transferring',0)
      num_working = 0
      for status in t_status.keys() :
        if ( status in ['idle','running','unsubmitted']) :
          num_working += t_status[status]
      if ( num_transferring < maxJob/2 and num_working < maxJob/2) :
        failed_jobs = self.collectFailed()
        for job in failed_jobs :
          if self.jobs_status.get(job).get('failed',0) == 0 : continue
          try :
            res = crabCommand('resubmit', job )
          except httplib.HTTPException, e:
            print "HTTPException",e," Another job will be resubmitted!"
            continue
          num_working += self.jobs_status.get(job).get('failed',0)
          if ( num_working> maxJob/2 ) : break
      time.sleep(3600)

      
if __name__ == "__main__":
  if ( len(sys.argv) !=2 ) :
    print "Wrong argument. Input just patten for jobs. egs) crab_v7-6-4"
    sys.exit(-1)
  cm = CrabMonitor(sys.argv[1])
  cm.Monitoring(20000)
  """
  cm.status()
  cm.printRemainJobsStatus()
  #cm.printCompletedJobs()
  totalStatus = cm.totalNumOfStatus()
  for status in totalStatus.keys() :
    print "Total %s : %d"%(status,totalStatus[status])
  """
