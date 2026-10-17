import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomTau3MuGunProducer2",
    PGunParameters = cms.PSet(
        MaxPtDs     = cms.double(120),
        MinPtDs     = cms.double(0.1),
        MaxEta     = cms.double(3.5),
        MinEta     = cms.double(1e-6),
        MaxPhi     = cms.double(3.141592653589793),
        MinPhi     = cms.double(-3.141592653589793),
        PartID     = cms.vint32(-13),
        RandomCharge = cms.bool(True)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('Tau3Mu decay'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
