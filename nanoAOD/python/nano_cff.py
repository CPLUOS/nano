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
    fileName = cms.untracked.string('nanoAOD.root')
    if hasattr(process, 'NANOAODSIMoutput'):
        # MC
        process.NANOAODSIMoutput.fileName = fileName
        jecFile = 'Summer16_07Aug2017_V10_MC'
    else:
        # DATA
        process.NANOAODoutput.fileName = fileName
        jecFile = 'Summer16_07Aug2017All_V10_DATA'
        
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

    # JEC
    from CondCore.CondDB.CondDB_cfi import CondDB
    if hasattr(CondDB, 'connect'): delattr(CondDB, 'connect')
    process.jec = cms.ESSource("PoolDBESSource",CondDB,
        connect = cms.string('sqlite_fip:nano/nanoAOD/data/jec/%s.db'%jecFile),            
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_%s_AK4PF"%jecFile),
                label= cms.untracked.string("AK4PF")),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFchs"%jecFile),
                label= cms.untracked.string("AK4PFchs")),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFPuppi"%jecFile),
                label= cms.untracked.string("AK4PFPuppi")),
            )
        )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
    print "JEC based on", process.jec.connect
        
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
    process.MessageLogger.cerr.FwkSummary.reportEvery = cms.untracked.int32(1000)
    
    return(process)
