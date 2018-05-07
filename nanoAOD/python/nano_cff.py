import FWCore.ParameterSet.Config as cms

def customiseMuons(process):
    # additional variables needed for h2mu
    process.muonTable.variables.globalMu = process.muonTable.variables.isPFcand
    process.muonTable.variables.globalMu.expr = cms.string('isGlobalMuon')
    process.muonTable.variables.globalMu.doc = process.muonTable.variables.globalMu.expr

    process.muonTable.variables.trackerMu = process.muonTable.variables.isPFcand
    process.muonTable.variables.trackerMu.expr = cms.string('isTrackerMuon')
    process.muonTable.variables.trackerMu.doc = process.muonTable.variables.trackerMu.expr

    return(process)
    
def customise(process, doHadron=True, fastSim=False):
    customiseMuons(process)
    
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    if doHadron:
        process.load('nano.nanoAOD.hadrons_cff')
        process.nanoAOD_step += process.hadTables

        process.load('nano.nanoAOD.v0_cff')
        process.nanoAOD_step += process.v0GenParticles + process.v0Tables
        
        from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
        run2_miniAOD_80XLegacy.toReplaceWith(process.nanoAOD_step, process.nanoAOD_step.copyAndExclude([process.v0GenParticles,process.v0Tables]))

    if fastSim:
        process.nanoAOD_step.remove(process.triggerObjectTable)
        process.nanoAOD_step.remove(process.l1bits)
    
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.MessageLogger.cerr.FwkSummary.reportEvery = cms.untracked.int32(1000)

    fileName = cms.untracked.string('nanoAOD.root')
    if hasattr(process, 'NANOAODSIMoutput'):          
        process.NANOAODSIMoutput.fileName = fileName
    else:
        process.NANOAODoutput.fileName = fileName
    
    return(process)
