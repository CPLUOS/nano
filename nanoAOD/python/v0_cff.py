import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

kshortCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedKshortVertices"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("Kshort"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
kshortCandidateTable.variables.pt.precision=14
kshortCandidateTable.variables.phi.precision=14
kshortCandidateTable.variables.eta.precision=14
kshortCandidateTable.variables.mass.precision=14

lambdaCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedLambdaVertices"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("Lambda"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
lambdaCandidateTable.variables.pt.precision=14
lambdaCandidateTable.variables.phi.precision=14
lambdaCandidateTable.variables.eta.precision=14
lambdaCandidateTable.variables.mass.precision=14

# Gen particle corresponding to Ks/Lambdaimport FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### User floats producers, selectors ##########################

v0GenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
	"drop *",
        "keep+ abs(pdgId) == 310 ",  #  keep first gen decay product for all kshort
        "keep+ abs(pdgId) == 3122 ",  #  keep first gen decay product for all lambda 		
   )
)

##################### Tables for final output and docs ##########################
v0GenParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("v0GenParticles"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("V0GenPart"),
    doc = cms.string("V0 gen particles "),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the taus
    variables = cms.PSet(
         pt  = Var("pt",  float,precision=8),
         phi = Var("phi", float,precision=8),
         eta  = Var("eta",  float,precision=8),
         mass = Var("mass", float,precision=8),
         pdgId  = Var("pdgId", int, doc="PDG id"),
         # x = Var("daughterRef(0).vx()", float, doc = "secondary vertex X position, in cm",precision=10),
         # y = Var("daughterRef(0).vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
         # z = Var("daughterRef(0).vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
    )
)

#before cross linking
v0Sequence = cms.Sequence(v0GenParticles)
#after cross linkining
v0Tables = cms.Sequence(kshortCandidateTable+lambdaCandidateTable+v0GenParticleTable)
