import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
jetMLIDTable = cms.EDProducer("JetMLIDProducer",
  jetLabel = cms.InputTag("linkedObjects","jets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  useQualityCuts = cms.bool(False),
)

jetMLIDTables = cms.Sequence(jetMLIDTable)

