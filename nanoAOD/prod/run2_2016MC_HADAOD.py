import sys

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: run2_2016MC -s NANO -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein file:ttbar_mc.root --conditions auto:run2_mc --era Run2_2016,run2_miniAOD_80XLegacy --customise nano/nanoAOD/nano_cff.customise
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2016,eras.run2_miniAOD_80XLegacy)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#process.maxEvents.input = cms.untracked.int32(1000)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/AODSIM/1001.root"),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('run2_2016MC nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('nanoAOD_AOD.root'),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#get all genparticles
process.load('PhysicsTools.NanoAOD.genparticles_cff')
from  PhysicsTools.NanoAOD.common_cff import *
process.genParticleTable.src = cms.InputTag("genParticles")
process.genParticleTable.variables.mass = Var("mass", float,precision=8,doc="Mass")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
process.load('nano.nanoAOD.hadrons_cff')
process.hadTable.jetLabel = cms.InputTag("patJets")
process.hadTable.vertexLabel = cms.InputTag("offlinePrimaryVertices")
process.hadTable.pfCandLabel = cms.InputTag("particleFlow")

process.pvTable =  cms.EDProducer("VertexTableProducer",
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    svSrc = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    svCut = cms.string(""),
    dlenMin = cms.double(0),
    dlenSigMin = cms.double(3),
    pvName = cms.string("PV"),
    svName = cms.string("SV"),
    svDoc  = cms.string("secondary vertices from IVF algorithm"),
)

process.load("Validation.RecoTrack.TrackValidation_cff")
#process.load('SimTracker.TrackerHitAssociation.tpClusterProducer_cfi')
#process.load('SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi')
#process.load('SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi')
#process.load('nano.nanoAOD.hadTruth_cff')

process.p = cms.Path(process.makePatJets+process.hadTables
                         #+process.genParticleTable
                         +process.mix+process.tracksValidationTruth
                         #+process.tpClusterProducer+process.quickTrackAssociatorByHits
                         #+process.trackingParticleRecoTrackAsssociation
#                         +process.hadTruthTables
	                 +process.pvTable
		    )

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.p,process.endjob_step,process.NANOAODSIMoutput_step)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
#from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC
#process = nanoAOD_customizeMC(process)

