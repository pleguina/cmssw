import FWCore.ParameterSet.Config as cms

# Reference:
#     https://github.com/cms-sw/genproductions/blob/master/genfragments/EightTeV/SingleMuMinusFlatPt0p2To100_cff.py

generator = cms.EDProducer("FlatRandomBeamHaloGunProducer2",
    PGunParameters = cms.PSet(
        MaxP = cms.double(120.0),
        MinP = cms.double(2.0),
        MaxVR = cms.double(700.0),
        MinVR = cms.double(200.0),
        MaxVZ = cms.double(1200.0),
        MinVZ = cms.double(1100.0),
        MaxTZ = cms.double(850.0),
        MinTZ = cms.double(750.0),
        MaxEta = cms.double(2.4),
        MinEta = cms.double(1.2),
        MaxPhi = cms.double(3.141592653589793),
        MinPhi = cms.double(-3.141592653589793),
        PartID = cms.vint32(-13),
        PSpectrum = cms.string('flatP'),
        RandomCharge = cms.bool(True)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single muon+/- p 2 to 120'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
