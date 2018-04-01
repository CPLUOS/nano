import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
hadTable = cms.EDProducer("HadronProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  pfCandLabel = cms.InputTag("packedPFCandidates"),
  applySoftLeptonCut = cms.bool(True),
  doFullMCMatching = cms.bool(False),
  # -- cuts on initial track collection --
  # Track normalized Chi2 <
  tkChi2Cut = cms.double(100),
  # Number of valid hits on track >=
  tkNHitsCut = cms.int32(0),
  # Pt of track >
  tkPtCut = cms.double(0.),
  # Track impact parameter significance >
  tkIPSigXYCut = cms.double(100000),
  tkIPSigZCut = cms.double(100000),
  
  # -- cuts on the vertex --
  # Vertex chi2 <
  vtxChi2Cut = cms.double(10),
  # XY decay distance significance >
  vtxDecaySigXYCut = cms.double(-1),
  # XYZ decay distance significance >
  vtxDecaySigXYZCut = cms.double(-1.),
  
  # -- miscellaneous cuts --
  # POCA distance between tracks <
  tkDCACut = cms.double(100),
  # cos(angleXY) between x and p of V0 candidate >
  cosThetaXYCut = cms.double(100),
  # cos(angleXYZ) between x and p of V0 candidate >
  cosThetaXYZCut = cms.double(100),
)

hadCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("hadTable",""),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("had"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)

hadCandidateTable.variables.pt.precision=14
hadCandidateTable.variables.phi.precision=14
hadCandidateTable.variables.eta.precision=14
hadCandidateTable.variables.mass.precision=14

hadJetTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("hadTable","jet"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("had_jet"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars,
        btagCMVA = Var("bDiscriminator('pfCombinedMVAV2BJetTags')",float,doc="CMVA V2 btag discriminator",precision=10),
        btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
        btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
        btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",float,doc="DeepCSV charm btag discriminator",precision=10),
        #puIdDisc = Var("userFloat('pileupJetId:fullDiscriminant')",float,doc="Pilup ID discriminant",precision=10),
        #puId = Var("userInt('pileupJetId:fullId')",int,doc="Pilup ID flags"),
        #jetId = Var("userInt('tightId')*2+userInt('looseId')",int,doc="Jet ID flags bit1 is loose, bit2 is tight"),
        #qgl = Var("userFloat('QGTagger:qgLikelihood')",float,doc="Quark vs Gluon likelihood discriminator",precision=10),
    ),
)

hadTables = cms.Sequence(hadTable+hadCandidateTable+hadJetTable)

