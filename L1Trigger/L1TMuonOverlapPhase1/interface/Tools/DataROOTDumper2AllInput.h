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
      const std::vector<int>& simTrackIdToGenIdx,   // simTrackId → genMuon 1-indexed (0=not a gen muon)
      unsigned int processorPhiZero,
      const std::vector<uint32_t>& stub_detId,
      const std::vector<short>& stub_phiHw,
      std::vector<signed char>& stub_trackId,
      std::vector<uint8_t>& stub_ambiguous);

  TTree* allInputTree = nullptr;

  // --- branches of OMTFAllInputTree ---
  unsigned int   reg_eventNum   = 0;
  unsigned char  reg_iProcessor = 0;
  signed char    reg_mtfType    = 0;  // l1t::tftype cast to int8 (omtf_pos or omtf_neg)

  std::vector<signed char> reg_stub_layer;    // logic layer index [0, nLayers)
  std::vector<short>       reg_stub_phiHw;    // absolute phi in HW units
  std::vector<short>       reg_stub_phiBHw;   // DT bending angle HW (0 for non-DT)
  std::vector<signed char> reg_stub_etaHw;    // eta HW (clamped to [-127, 127])
  std::vector<short>       reg_stub_r;        // radial distance [cm]
  std::vector<signed char> reg_stub_quality;  // qualityHw
  std::vector<signed char> reg_stub_type;     // MuonStub::Type cast to int8
  std::vector<signed char> reg_stub_bx;       // BX offset (clamped to [-127, 127])
  // SimTrack truth (only filled when doSimTruth_ == true)
  std::vector<signed char> reg_stub_trackId;   // 0 = noise, 1..K = gen-muon idx
  std::vector<uint8_t>     reg_stub_ambiguous; // 1 if dominant track < 50% digis

  // --- per-event cache ---
  // Filled per-processor in observeProcesorEmulation, flushed in observeEventEnd.
  struct CachedRegion {
    unsigned char iProcessor = 0;
    signed char   mtfType    = 0;
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
  edm::InputTag simTrackTag_;
  edm::InputTag genParticleTag_;
};

#endif /* L1T_OmtfP1_TOOLS_DATAROOTDUMPER2ALLINPUT_H_ */
