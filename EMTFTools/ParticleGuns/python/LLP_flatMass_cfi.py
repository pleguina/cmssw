import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomLLPGunProducer2",
    PGunParameters = cms.PSet(
        MaxMassH   = cms.double(1000),
        MinMassH   = cms.double(20),
        MaxPtH     = cms.double(120),
        MinPtH     = cms.double(1),
        MaxCTauLLP = cms.double(5000),
        MinCTauLLP = cms.double(10),
        MaxEta     = cms.double(3.5),
        MinEta     = cms.double(1e-6),
        MaxPhi     = cms.double(3.141592653589793),
        MinPhi     = cms.double(-3.141592653589793),
        LLPMassSpectrum = cms.string('flatMass'),
        PartID     = cms.vint32(-13)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('LLP decay with flat mass spectrum'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
