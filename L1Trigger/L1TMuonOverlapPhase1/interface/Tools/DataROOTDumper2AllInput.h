/*
 * DataROOTDumper2AllInput.h
 *
 * Extension of DataROOTDumper2 that writes a second ROOT tree (OMTFAllInputTree)
 * containing ALL input stubs for each processor region per event — not only the
 * stubs matched by the best golden pattern.
 *
 * This is needed to train the Object Condensation and noise-rejection heads of the
 * HECIN+OC GNN, which require the full view of the detector hits in a region.
 *
 * The existing OMTFHitsTree (per-candidate stubs, from the base class) is still
 * written unchanged. Join the two trees on (reg_eventNum, reg_iProcessor,
 * reg_mtfType) to correlate regions with candidates.
 *
 * Author: pleguina (extended from DataROOTDumper2 by kbunkow)
 */

#ifndef L1T_OmtfP1_TOOLS_DATAROOTDUMPER2ALLINPUT_H_
#define L1T_OmtfP1_TOOLS_DATAROOTDUMPER2ALLINPUT_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper2.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EmulationObserverBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/AngleConverterBase.h"  // MuonGeometryTokens

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include <vector>
#include <memory>

class TTree;

// ---------------------------------------------------------------------------
// DataROOTDumper2AllInput
// ---------------------------------------------------------------------------
// Adds a second ROOT tree (OMTFAllInputTree) on top of the existing OMTFHitsTree.
// Each entry of OMTFAllInputTree corresponds to one (event, processor-region) pair
// and holds ALL input stubs to that region (before pattern matching).
//
// Branch naming convention: reg_* for per-region scalars, reg_stub_* for the
// parallel stub vectors (same index → same stub).
//
// Two SimTrack-level truth branches are added when the digi-sim link
// InputTags are configured (i.e. doSimTruth == true):
//   reg_stub_trackId   : 0 = noise/unmatched, 1..K = gen-muon index (1-indexed)
//   reg_stub_ambiguous : 1 if the dominant SimTrack covers < 50 % of matched digis
// ---------------------------------------------------------------------------
class DataROOTDumper2AllInput : public DataROOTDumper2 {
public:
  DataROOTDumper2AllInput(const edm::ParameterSet& edmCfg,
                          const OMTFConfiguration* omtfConfig,
                          CandidateSimMuonMatcher* candidateSimMuonMatcher,
                          const MuonGeometryTokens& muonGeometryTokens);

  ~DataROOTDumper2AllInput() override;

  void beginRun(edm::EventSetup const& eventSetup) override;

  // Overrides the no-op in DataROOTDumper2: loops all layers/stubs in the
  // OMTFinput and caches them for writing in observeEventEnd.
  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>& input,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const FinalMuons& finalMuons) override;

  // Flushes the per-event region cache to OMTFAllInputTree, then calls the
  // base-class routine which writes OMTFHitsTree as before.
  void observeEventEnd(const edm::Event& iEvent,
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override;

  void endJob() override;

private:
  void initializeAllInputTree();

  // Assign reg_stub_trackId / reg_stub_ambiguous for one cached region.
  // Returns false (silently) when digi-sim link handles are not available.
  void assignSimTruth(
      const edm::Event& iEvent,
      const edm::Handle<MuonDigiCollection<DTLayerId, DTDigiSimLink>>& dtLinksH,
      const edm::Handle<edm::DetSetVector<RPCDigiSimLink>>& rpcLinksH,
      const edm::Handle<edm::DetSetVector<StripDigiSimLink>>& cscLinksH,
      const edm::Handle<edm::PSimHitContainer>& dtSimHitsH,
      const edm::Handle<edm::PSimHitContainer>& rpcSimHitsH,
      const edm::Handle<edm::PSimHitContainer>& cscSimHitsH,
      const std::vector<int>& simTrackIdToGenIdx,   // simTrackId → genMuon 1-indexed (0=not a gen muon)
      const std::vector<int>& simTrackIdToTpIdx,    // simTrackId → TpMuon 1-indexed (0=not a selected TP muon)
      unsigned int processorPhiZero,
      const std::vector<uint32_t>& stub_detId,
      const std::vector<short>& stub_phiHw,
      std::vector<signed char>& stub_trackId,
      std::vector<short>& stub_tpId,
      std::vector<uint8_t>& stub_ambiguous,
      std::vector<float>& stub_tof,
      std::vector<float>& stub_tofSpread,
      std::vector<signed char>& stub_nSimHit);

  TTree* allInputTree = nullptr;

  // --- branches of OMTFAllInputTree ---
  unsigned int   reg_runNum     = 0;
  unsigned int   reg_lumiNum    = 0;
  unsigned long long reg_eventNum64 = 0;
  unsigned int   reg_eventNum   = 0;
  unsigned char  reg_iProcessor = 0;
  signed char    reg_mtfType    = 0;  // l1t::tftype cast to int8 (omtf_pos or omtf_neg)
  signed char    reg_endcap     = 0;  // +1 = OMTF_POS, -1 = OMTF_NEG, 0 = unknown

  std::vector<signed char> reg_stub_layer;    // logic layer index [0, nLayers)
  std::vector<short>       reg_stub_phiHw;    // absolute phi in HW units
  std::vector<short>       reg_stub_phiBHw;   // DT bending angle HW (0 for non-DT)
  std::vector<signed char> reg_stub_etaHw;    // eta HW (clamped to [-127, 127])
  std::vector<short>       reg_stub_r;        // radial distance [cm]
  std::vector<signed char> reg_stub_quality;  // qualityHw
  std::vector<signed char> reg_stub_type;     // MuonStub::Type cast to int8
  std::vector<signed char> reg_stub_bx;       // BX offset (clamped to [-127, 127])
  std::vector<uint32_t>    reg_stub_detId;    // stable source chamber/raw detId per stub
  // SimTrack truth (only filled when doSimTruth_ == true)
  std::vector<signed char> reg_stub_trackId;   // 0 = noise, 1..K = gen-muon idx
  std::vector<short>       reg_stub_tpId;      // 0 = no selected TP muon, 1..M = TpMuon idx
  std::vector<uint8_t>     reg_stub_ambiguous; // 1 if dominant track < 50% digis
  std::vector<float>       reg_stub_tof;       // mean SimHit TOF [ns], -999 when unavailable
  std::vector<float>       reg_stub_tofSpread; // max-min SimHit TOF [ns], 0 when unavailable
  std::vector<signed char> reg_stub_nSimHit;   // number of SimHits used for TOF (clamped to 127)

  // Selected TrackingParticle muons persisted per event (duplicated for each region entry).
  std::vector<float>       TpMuon_pt;
  std::vector<float>       TpMuon_eta;
  std::vector<float>       TpMuon_phi;
  std::vector<signed char> TpMuon_charge;
  std::vector<float>       TpMuon_dxy;
  std::vector<float>       TpMuon_lxy;
  // Raw production-vertex transverse coordinates [cm]. Persisted alongside the
  // standard straight-line TpMuon_dxy (TrackingParticle::dxy()) so that a
  // curvature-consistent (helix-based) impact parameter can be reconstructed
  // downstream using TpMuon_pt/TpMuon_charge/TpMuon_phi and the known B field.
  // See DATASETS_INFO.txt "KNOWN CAVEAT: TpMuon_dxy is NOT the gun's true d0".
  std::vector<float>       TpMuon_vx;
  std::vector<float>       TpMuon_vy;
  // Curvature-consistent (helix-based) impact parameter, computed the same way
  // as CandidateSimMuonMatcher::MatchingResult's SimTrack-path constructor
  // (rg from qinvpt/B, then d0 = rg + q*|center|), but applied to the
  // TrackingParticle-path muons, which that class still leaves on the naive
  // straight-line trackingParticle.dxy(). This is the value that should
  // actually stay within [MinDxy, MaxDxy] of the gun.
  std::vector<float>       TpMuon_dxyHelix;
  std::vector<int>         TpMuon_eventId;
  std::vector<signed char> TpMuon_bunchCrossing;
  std::vector<uint8_t>     TpMuon_isPileup;
  std::vector<uint8_t>     TpMuon_isInTime;
  std::vector<signed char> TpMuon_originClass;  // 0=HS, 1=PU-in-time, 2=PU-out-of-time, 3=other

  int reg_pu_nTpMuonInTime = 0;
  int reg_pu_nTpMuonOOT = 0;

  // --- per-event cache ---
  // Filled per-processor in observeProcesorEmulation, flushed in observeEventEnd.
  struct CachedRegion {
    unsigned char iProcessor = 0;
    signed char   mtfType    = 0;
    signed char   endcap     = 0;  // +1 = OMTF_POS, -1 = OMTF_NEG, 0 = unknown
    std::vector<signed char> stub_layer;
    std::vector<short>       stub_phiHw;
    std::vector<short>       stub_phiBHw;
    std::vector<signed char> stub_etaHw;
    std::vector<short>       stub_r;
    std::vector<signed char> stub_quality;
    std::vector<signed char> stub_type;
    std::vector<signed char> stub_bx;
    // cached for sim-truth assignment in observeEventEnd
    std::vector<uint32_t>    stub_detId;
  };

  std::vector<CachedRegion> cachedRegions;  // up to ~12 entries per event

  // --- SimTrack truth infrastructure ---
  bool doSimTruth_ = false;

  const MuonGeometryTokens& muonGeometryTokens_;
  edm::ESWatcher<MuonGeometryRecord> geomWatcher_;
  edm::ESHandle<RPCGeometry> georpc_;
  edm::ESHandle<CSCGeometry> geocsc_;
  edm::ESHandle<DTGeometry>  geodt_;

  edm::InputTag dtSimLinkTag_;
  edm::InputTag rpcSimLinkTag_;
  edm::InputTag cscSimLinkTag_;
  edm::InputTag dtSimHitTag_;
  edm::InputTag rpcSimHitTag_;
  edm::InputTag cscSimHitTag_;
  edm::InputTag simTrackTag_;
  edm::InputTag genParticleTag_;
  edm::InputTag trackingParticleTag_;
};

#endif /* L1T_OmtfP1_TOOLS_DATAROOTDUMPER2ALLINPUT_H_ */
