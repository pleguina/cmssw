"""
G3 — Single Displaced Muon, OMTF overlap region, no PU.

GMT-visible-stub campaign: clean 1-candidate displaced overlap reference.

Key properties:
  n muons/event: 1  (PartID has exactly ONE entry; charge randomised by gun)
  pT: flat in 1/pT over [2, 200] GeV -> unbiased curvature coverage
  eta: 0.82 < |eta| < 1.24  (OMTF overlap acceptance)
       generated in [-1.24, +1.24]; etaFilter enforces |eta| in [0.82, 1.24]
  phi: full 2pi
  displacement: flat d0 in [0, 50] cm; MaxLxy = 200 cm (vertex inside MB1)
  PU: none

Charge note:
  PartID = [-13] (mu-) with RandomCharge=True -> ~50% mu+ / ~50% mu-.
  Exactly one PDG ID entry satisfies the single-muon rule.

Eta filter note:
  etaFilter is defined AND included in ProductionFilterSequence.

Production target: 250,000 events / 500 jobs at 500 events/job.
"""
import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("FlatRandomPtGunProducer2",
    PGunParameters = cms.PSet(
        PartID         = cms.vint32(-13),            # single muon; charge randomised below
        MinPt          = cms.double(2.0),             # [GeV] pT range
        MaxPt          = cms.double(200.0),
        MinDxy         = cms.double(0.0),             # [cm] d0 range, flat uniform
        MaxDxy         = cms.double(50.0),            # [cm] full OMTF displaced coverage
        MaxLxy         = cms.double(200.0),           # [cm] keep vertex inside MB1 (r < 231 cm)
        MinEta         = cms.double(-1.24),           # full OMTF acceptance; filter narrows
        MaxEta         = cms.double( 1.24),
        MinPhi         = cms.double(-3.14159265359),
        MaxPhi         = cms.double( 3.14159265359),
        PtSpectrum     = cms.string('flatOneOverPt'),
        VertexSpectrum = cms.string('flatD0'),
        RandomCharge   = cms.bool(True),              # 50% mu+ / 50% mu- per event
    ),
    Verbosity       = cms.untracked.int32(0),
    psethack        = cms.string('single displaced muon flatOneOverPt flatD0 d0 0-50cm OMTF overlap'),
    AddAntiParticle = cms.bool(False),
    firstRun        = cms.untracked.uint32(1),
)

# Eta filter: keep only events where muon is strictly in OMTF overlap acceptance
etaFilter = cms.EDFilter("MCSingleParticleFilter",
    ParticleID  = cms.untracked.vint32(13, -13, 13, -13),
    Status      = cms.untracked.vint32(1, 1, 1, 1),
    MinPt       = cms.untracked.vdouble(0.0, 0.0, 0.0, 0.0),
    MinEta      = cms.untracked.vdouble(0.82, 0.82, -1.24, -1.24),
    MaxEta      = cms.untracked.vdouble(1.24, 1.24, -0.82, -0.82),
)

ProductionFilterSequence = cms.Sequence(generator * etaFilter)
