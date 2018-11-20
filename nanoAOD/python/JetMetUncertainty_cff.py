import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

jetUncEvaluatorTable = cms.EDProducer("JetUncertaintyEvaluator",
  src = cms.InputTag("linkedObjects","jets"),
  rho = cms.InputTag('fixedGridRhoFastjetAll'), 
  payloadName  = cms.string('AK4PFchs'), 
  jetResFile   = cms.string("nano/nanoAOD/data/jer/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"), 
  jetResSFFile = cms.string("nano/nanoAOD/data/jer/Summer16_25nsV1_MC_SF_AK4PFchs.txt"), 
)

jetUncEvaluatorTables = cms.Sequence(jetUncEvaluatorTable)

