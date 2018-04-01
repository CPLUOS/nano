import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
hadTruthTable = cms.EDProducer("HadTruthProducer",
  recoRecoToSim = cms.InputTag("trackingParticleRecoTrackAsssociation"),
  recoSimToReco = cms.InputTag("trackingParticleRecoTrackAsssociation"),
  hadronCands = cms.InputTag("hadTable"),
  genLabel  = cms.InputTag("genParticles"),
  trackingVertexLabel = cms.InputTag("mix", "MergedTrackTruth"),
  trackingParticleLabel = cms.InputTag("mix", "MergedTrackTruth")
)

hadTruthCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("hadTruthTable"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("genHadron"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=14),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=14),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)

hadTruthCandidateTable.variables.pt.precision=14
hadTruthCandidateTable.variables.phi.precision=14
hadTruthCandidateTable.variables.eta.precision=14
hadTruthCandidateTable.variables.mass.precision=14

hadTruthTables = cms.Sequence(hadTruthTable+hadTruthCandidateTable)
