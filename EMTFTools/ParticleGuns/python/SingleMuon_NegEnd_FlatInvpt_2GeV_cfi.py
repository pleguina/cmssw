import FWCore.ParameterSet.Config as cms

# Reference:
#     https://github.com/cms-sw/genproductions/blob/master/genfragments/EightTeV/SingleMuMinusFlatPt0p2To100_cff.py

generator = cms.EDProducer("FlatRandomPtGunProducer2",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(120.0),
        MinPt = cms.double(2.0),
        MaxDxy = cms.double(120.0),
        MinDxy = cms.double(0.0),
        MaxLxy = cms.double(300.0),
        MaxEta = cms.double(-1.2),
        MinEta = cms.double(-2.5),
        MaxPhi = cms.double(3.141592653589793),
        MinPhi = cms.double(-3.141592653589793),
        PartID = cms.vint32(-13),
        PtSpectrum = cms.string('flatOneOverPt'),
        VertexSpectrum = cms.string('none'),
        RandomCharge = cms.bool(True)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single muon+/- pt 2 to 120 flat in 1/pt negative endcap'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
