/*
 * DataROOTDumper2AllInput.cc
 *
 * See DataROOTDumper2AllInput.h for description.
 *
 * Author: pleguina
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper2AllInput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TTree.h"

#include <cmath>
#include <limits>
#include <map>
#include <set>

// ---------------------------------------------------------------------------
DataROOTDumper2AllInput::DataROOTDumper2AllInput(const edm::ParameterSet& edmCfg,
                                                 const OMTFConfiguration* omtfConfig,
                                                 CandidateSimMuonMatcher* candidateSimMuonMatcher,
                                                 const MuonGeometryTokens& muonGeometryTokens)
    : DataROOTDumper2(edmCfg, omtfConfig, candidateSimMuonMatcher),
      muonGeometryTokens_(muonGeometryTokens) {

  // Sim-truth is optional: only enabled when the dtDigiSimLinksInputTag parameter exists
  if (edmCfg.exists("dtDigiSimLinksInputTag")) {
    doSimTruth_ = true;

    dtSimLinkTag_  = edmCfg.getParameter<edm::InputTag>("dtDigiSimLinksInputTag");
    rpcSimLinkTag_ = edmCfg.getParameter<edm::InputTag>("rpcDigiSimLinkInputTag");
    cscSimLinkTag_ = edmCfg.getParameter<edm::InputTag>("cscStripDigiSimLinksInputTag");
    dtSimHitTag_   = edmCfg.exists("dtSimHitsInputTag")
               ? edmCfg.getParameter<edm::InputTag>("dtSimHitsInputTag")
               : edm::InputTag("g4SimHits", "MuonDTHits");
    rpcSimHitTag_  = edmCfg.exists("rpcSimHitsInputTag")
               ? edmCfg.getParameter<edm::InputTag>("rpcSimHitsInputTag")
               : edm::InputTag("g4SimHits", "MuonRPCHits");
    cscSimHitTag_  = edmCfg.exists("cscSimHitsInputTag")
               ? edmCfg.getParameter<edm::InputTag>("cscSimHitsInputTag")
               : edm::InputTag("g4SimHits", "MuonCSCHits");
    simTrackTag_   = edmCfg.getParameter<edm::InputTag>("simTracksTag");
    genParticleTag_= edmCfg.getParameter<edm::InputTag>("genParticleTag");

    edm::LogVerbatim("l1tOmtfEventPrint") << "DataROOTDumper2AllInput: doSimTruth = true" << std::endl;
  }

  initializeAllInputTree();
  edm::LogVerbatim("l1tOmtfEventPrint") << "DataROOTDumper2AllInput created" << std::endl;
}

DataROOTDumper2AllInput::~DataROOTDumper2AllInput() {}

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::beginRun(edm::EventSetup const& eventSetup) {
  DataROOTDumper2::beginRun(eventSetup);

  if (doSimTruth_ && geomWatcher_.check(eventSetup)) {
    georpc_ = eventSetup.getHandle(muonGeometryTokens_.rpcGeometryEsToken);
    geocsc_ = eventSetup.getHandle(muonGeometryTokens_.cscGeometryEsToken);
    geodt_  = eventSetup.getHandle(muonGeometryTokens_.dtGeometryEsToken);
  }
}

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::initializeAllInputTree() {
  edm::Service<TFileService> fs;

  allInputTree =
      fs->make<TTree>("OMTFAllInputTree", "All input stubs per processor region per event");

  allInputTree->Branch("reg_eventNum",    &reg_eventNum,  "reg_eventNum/i");
  allInputTree->Branch("reg_iProcessor",  &reg_iProcessor, "reg_iProcessor/b");
  allInputTree->Branch("reg_mtfType",     &reg_mtfType, "reg_mtfType/B");

  allInputTree->Branch("reg_stub_layer",   &reg_stub_layer);
  allInputTree->Branch("reg_stub_phiHw",   &reg_stub_phiHw);
  allInputTree->Branch("reg_stub_phiBHw",  &reg_stub_phiBHw);
  allInputTree->Branch("reg_stub_etaHw",   &reg_stub_etaHw);
  allInputTree->Branch("reg_stub_r",       &reg_stub_r);
  allInputTree->Branch("reg_stub_quality", &reg_stub_quality);
  allInputTree->Branch("reg_stub_type",    &reg_stub_type);
  allInputTree->Branch("reg_stub_bx",      &reg_stub_bx);

  // SimTrack truth branches are always created so the tree schema is stable.
  // When doSimTruth_ == false they stay filled with zeros.
  allInputTree->Branch("reg_stub_trackId",   &reg_stub_trackId);
  allInputTree->Branch("reg_stub_ambiguous", &reg_stub_ambiguous);
  allInputTree->Branch("reg_stub_tof",       &reg_stub_tof);
  allInputTree->Branch("reg_stub_tofSpread", &reg_stub_tofSpread);
  allInputTree->Branch("reg_stub_nSimHit",   &reg_stub_nSimHit);
}

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::observeProcesorEmulation(unsigned int iProcessor,
                                                       l1t::tftype mtfType,
                                                       const std::shared_ptr<OMTFinput>& omtfInput,
                                                       const AlgoMuons& /*algoCandidates*/,
                                                       const AlgoMuons& /*gbCandidates*/,
                                                       const FinalMuons& /*finalMuons*/) {
  CachedRegion region;
  region.iProcessor = static_cast<unsigned char>(iProcessor);
  region.mtfType    = static_cast<signed char>(mtfType);

  // getMuonStubs() returns MuonStubPtrs2D = vector<vector<MuonStubPtr>>
  // indexed [iLayer][iStub]. Iterate all layers and all stubs.
  const auto& muonStubsInLayers = omtfInput->getMuonStubs();
  for (unsigned int iLayer = 0; iLayer < muonStubsInLayers.size(); ++iLayer) {
    for (const auto& stub : muonStubsInLayers[iLayer]) {
      if (!stub)
        continue;

      region.stub_layer.push_back(static_cast<signed char>(iLayer));
      region.stub_phiHw.push_back(static_cast<short>(stub->phiHw));
      region.stub_phiBHw.push_back(static_cast<short>(stub->phiBHw));

      int etaHw = stub->etaHw;
      if (etaHw > 127)
        etaHw = 127;
      else if (etaHw < -127)
        etaHw = -127;
      region.stub_etaHw.push_back(static_cast<signed char>(etaHw));

      region.stub_r.push_back(static_cast<short>(stub->r));
      region.stub_quality.push_back(static_cast<signed char>(stub->qualityHw));
      region.stub_type.push_back(static_cast<signed char>(stub->type));

      int bx = stub->bx;
      if (bx > 127)
        bx = 127;
      else if (bx < -127)
        bx = -127;
      region.stub_bx.push_back(static_cast<signed char>(bx));

      // Store raw DetId for sim-truth lookup in observeEventEnd
      region.stub_detId.push_back(static_cast<uint32_t>(stub->detId));
    }
  }

  // Only store regions that actually had at least one stub
  if (!region.stub_layer.empty())
    cachedRegions.push_back(std::move(region));
}

// ---------------------------------------------------------------------------
// Helper: match one stub to SimTracks via digi-sim links.
// Returns the dominant SimTrack ID (or 0) and sets ambiguous.
// ---------------------------------------------------------------------------
namespace {

  // Build a map  simTrackId → count of matched digis  for one stub
  // using the three digi-sim link collections.
  void accumRPC(uint32_t detIdRaw,
                double stubGlobalPhi,
                const edm::DetSetVector<RPCDigiSimLink>& rpcLinks,
                const RPCGeometry& georpc,
                std::map<uint32_t, unsigned>& votes) {
    RPCDetId rpcDetId(detIdRaw);
    const RPCRoll* roll = georpc.roll(rpcDetId);
    if (!roll)
      return;
    const int nstripsRpc = roll->nstrips();
    // Use equal_range instead of find() to avoid assertion on duplicate DetSets
    auto rpcRange = std::equal_range(rpcLinks.begin(), rpcLinks.end(), edm::det_id_type(detIdRaw));
    for (auto dsIt = rpcRange.first; dsIt != rpcRange.second; ++dsIt) {
      for (const auto& link : *dsIt) {
        int strip = static_cast<int>(link.getStrip());
        if (strip < 1 || strip > nstripsRpc)
          continue;
        double phi = (roll->toGlobal(roll->centreOfStrip(strip))).phi();
        if (std::abs(stubGlobalPhi - phi) < 0.02)
          votes[link.getTrackId()]++;
      }
    }
  }

  void accumDT(uint32_t detIdRaw,
               double stubGlobalPhi,
               const MuonDigiCollection<DTLayerId, DTDigiSimLink>& dtLinks,
               const DTGeometry& geodt,
               std::map<uint32_t, unsigned>& votes) {
    const DTChamber* chamber = geodt.chamber(DTLayerId(detIdRaw));
    if (!chamber)
      return;
    for (auto* superlayer : chamber->superLayers()) {
      if (superlayer->id().superLayer() == 2)  // skip theta
        continue;
      for (auto* layer : superlayer->layers()) {
        const int firstWire = layer->specificTopology().firstChannel();
        const int lastWire  = layer->specificTopology().lastChannel();
        auto range = dtLinks.get(layer->id());
        for (auto lnk = range.first; lnk != range.second; ++lnk) {
          int wire = static_cast<int>(lnk->wire());
          if (wire < firstWire || wire > lastWire)
            continue;
          auto wireX   = layer->specificTopology().wirePosition(wire);
          auto digiPhi = layer->toGlobal(LocalPoint(wireX, 0, 0)).phi();
          if (std::abs(stubGlobalPhi - digiPhi) < 0.03)
            votes[lnk->SimTrackId()]++;
        }
      }
    }
  }

  void accumCSC(uint32_t detIdRaw,
                double stubGlobalPhi,
                const edm::DetSetVector<StripDigiSimLink>& cscLinks,
                const CSCGeometry& geocsc,
                std::map<uint32_t, unsigned>& votes) {
    const CSCChamber* chamber = geocsc.chamber(CSCDetId(detIdRaw));
    if (!chamber)
      return;
    for (const auto* layer : chamber->layers()) {
      // Use equal_range instead of find() to avoid assertion on duplicate DetSets
      auto cscRange = std::equal_range(
          cscLinks.begin(), cscLinks.end(), edm::det_id_type(layer->id().rawId()));
      const int nstrips = layer->geometry()->numberOfStrips();
      for (auto dsIt = cscRange.first; dsIt != cscRange.second; ++dsIt) {
        for (const auto& link : *dsIt) {
          int strip = static_cast<int>(link.channel());
          if (strip < 1 || strip > nstrips)
            continue;
          auto phi = layer->centerOfStrip(strip).phi();
          if (std::abs(stubGlobalPhi - phi) < 0.03)
            votes[link.SimTrackId()]++;
        }
      }
    }
  }

  void collectRPCTof(uint32_t detIdRaw,
                     double stubGlobalPhi,
                     const edm::PSimHitContainer& rpcSimHits,
                     const RPCGeometry& georpc,
                     std::vector<float>& tofs) {
    RPCDetId rpcDetId(detIdRaw);
    const RPCRoll* roll = georpc.roll(rpcDetId);
    if (!roll)
      return;
    for (const auto& simHit : rpcSimHits) {
      if (simHit.detUnitId() != detIdRaw)
        continue;
      const int strip = roll->strip(simHit.localPosition());
      if (strip < 1 || strip > roll->nstrips())
        continue;
      const double phi = roll->toGlobal(roll->centreOfStrip(strip)).phi();
      if (std::abs(stubGlobalPhi - phi) < 0.02)
        tofs.push_back(simHit.timeOfFlight());
    }
  }

  void collectDTTof(uint32_t detIdRaw,
                    double stubGlobalPhi,
                    const edm::PSimHitContainer& dtSimHits,
                    const DTGeometry& geodt,
                    std::vector<float>& tofs) {
    const DTChamber* chamber = geodt.chamber(DTLayerId(detIdRaw));
    if (!chamber)
      return;
    for (const auto& simHit : dtSimHits) {
      const DTLayer* layer = geodt.layer(DTLayerId(simHit.detUnitId()));
      if (!layer)
        continue;
      if (layer->chamber()->id().rawId() != detIdRaw)
        continue;
      const auto globalPoint = layer->toGlobal(simHit.localPosition());
      if (std::abs(stubGlobalPhi - globalPoint.phi()) < 0.03)
        tofs.push_back(simHit.timeOfFlight());
    }
  }

  void collectCSCTof(uint32_t detIdRaw,
                     double stubGlobalPhi,
                     const edm::PSimHitContainer& cscSimHits,
                     const CSCGeometry& geocsc,
                     std::vector<float>& tofs) {
    const CSCChamber* chamber = geocsc.chamber(CSCDetId(detIdRaw));
    if (!chamber)
      return;
    for (const auto& simHit : cscSimHits) {
      const CSCLayer* layer = geocsc.layer(CSCDetId(simHit.detUnitId()));
      if (!layer)
        continue;
      if (layer->chamber()->id().rawId() != detIdRaw)
        continue;
      const double strip = layer->geometry()->strip(simHit.localPosition());
      const double phi = layer->centerOfStrip(std::round(strip)).phi().value();
      if (std::abs(stubGlobalPhi - phi) < 0.03)
        tofs.push_back(simHit.timeOfFlight());
    }
  }

  void summarizeTof(const std::vector<float>& tofs, float& meanTof, float& spreadTof, signed char& nTof) {
    if (tofs.empty()) {
      meanTof = -999.f;
      spreadTof = 0.f;
      nTof = 0;
      return;
    }

    float sum = 0.f;
    float minTof = std::numeric_limits<float>::max();
    float maxTof = std::numeric_limits<float>::lowest();
    for (float t : tofs) {
      sum += t;
      if (t < minTof) minTof = t;
      if (t > maxTof) maxTof = t;
    }
    meanTof = sum / static_cast<float>(tofs.size());
    spreadTof = maxTof - minTof;
    const size_t capped = std::min<size_t>(tofs.size(), 127);
    nTof = static_cast<signed char>(capped);
  }

}  // anonymous namespace

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::assignSimTruth(
    const edm::Event& /*iEvent*/,
    const edm::Handle<MuonDigiCollection<DTLayerId, DTDigiSimLink>>& dtLinksH,
    const edm::Handle<edm::DetSetVector<RPCDigiSimLink>>& rpcLinksH,
    const edm::Handle<edm::DetSetVector<StripDigiSimLink>>& cscLinksH,
  const edm::Handle<edm::PSimHitContainer>& dtSimHitsH,
  const edm::Handle<edm::PSimHitContainer>& rpcSimHitsH,
  const edm::Handle<edm::PSimHitContainer>& cscSimHitsH,
    const std::vector<int>& simTrackIdToGenIdx,
    unsigned int processorPhiZero,
    const std::vector<uint32_t>& stub_detId,
    const std::vector<short>& stub_phiHw,
    std::vector<signed char>& stub_trackId,
  std::vector<uint8_t>& stub_ambiguous,
  std::vector<float>& stub_tof,
  std::vector<float>& stub_tofSpread,
  std::vector<signed char>& stub_nSimHit) {

  const unsigned int nStubs = stub_detId.size();
  stub_trackId.assign(nStubs, 0);
  stub_ambiguous.assign(nStubs, 0);
  stub_tof.assign(nStubs, -999.f);
  stub_tofSpread.assign(nStubs, 0.f);
  stub_nSimHit.assign(nStubs, 0);

  for (unsigned int iStub = 0; iStub < nStubs; ++iStub) {
    uint32_t detIdRaw = stub_detId[iStub];
    DetId detId(detIdRaw);
    if (detId.det() != DetId::Muon)
      continue;

    double globalPhi = omtfConfig->procHwPhiToGlobalPhi(stub_phiHw[iStub],
                                                         static_cast<int>(processorPhiZero));

    std::map<uint32_t, unsigned> votes;

    switch (detId.subdetId()) {
      case MuonSubdetId::RPC:
        accumRPC(detIdRaw, globalPhi, *rpcLinksH, *georpc_, votes);
        if (rpcSimHitsH.isValid()) {
          std::vector<float> tofs;
          collectRPCTof(detIdRaw, globalPhi, *rpcSimHitsH, *georpc_, tofs);
          summarizeTof(tofs, stub_tof[iStub], stub_tofSpread[iStub], stub_nSimHit[iStub]);
        }
        break;
      case MuonSubdetId::DT:
        accumDT(detIdRaw, globalPhi, *dtLinksH, *geodt_, votes);
        if (dtSimHitsH.isValid()) {
          std::vector<float> tofs;
          collectDTTof(detIdRaw, globalPhi, *dtSimHitsH, *geodt_, tofs);
          summarizeTof(tofs, stub_tof[iStub], stub_tofSpread[iStub], stub_nSimHit[iStub]);
        }
        break;
      case MuonSubdetId::CSC:
        accumCSC(detIdRaw, globalPhi, *cscLinksH, *geocsc_, votes);
        if (cscSimHitsH.isValid()) {
          std::vector<float> tofs;
          collectCSCTof(detIdRaw, globalPhi, *cscSimHitsH, *geocsc_, tofs);
          summarizeTof(tofs, stub_tof[iStub], stub_tofSpread[iStub], stub_nSimHit[iStub]);
        }
        break;
      default:
        break;
    }

    if (votes.empty())
      continue;

    // Find dominant SimTrack
    uint32_t bestId  = 0;
    unsigned bestCnt = 0;
    unsigned total   = 0;
    for (const auto& kv : votes) {
      total += kv.second;
      if (kv.second > bestCnt) {
        bestCnt = kv.second;
        bestId  = kv.first;
      }
    }

    // Map SimTrack ID to gen-muon 1-indexed (0 = not a gen muon in acceptance)
    int genIdx = 0;
    if (bestId < simTrackIdToGenIdx.size())
      genIdx = simTrackIdToGenIdx[bestId];

    stub_trackId[iStub]   = static_cast<signed char>(genIdx);
    // ambiguous if best track covers less than half the matched digis
    stub_ambiguous[iStub] = (total > 0 && bestCnt * 2 < total) ? 1 : 0;
  }
}

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::observeEventEnd(
    const edm::Event& iEvent,
    std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {

  const unsigned int eventNum = iEvent.id().event();

  // --- Optionally build SimTrack → gen-muon index map ---
  // simTrackIdToGenIdx[trackId] = 1-indexed gen-muon index, or 0 if not a gen muon
  std::vector<int> simTrackIdToGenIdx;   // indexed by SimTrack::trackId()

  edm::Handle<MuonDigiCollection<DTLayerId, DTDigiSimLink>> dtLinksH;
  edm::Handle<edm::DetSetVector<RPCDigiSimLink>>            rpcLinksH;
  edm::Handle<edm::DetSetVector<StripDigiSimLink>>          cscLinksH;
  edm::Handle<edm::PSimHitContainer>                        dtSimHitsH;
  edm::Handle<edm::PSimHitContainer>                        rpcSimHitsH;
  edm::Handle<edm::PSimHitContainer>                        cscSimHitsH;
  bool simHandlesOk = false;

  if (doSimTruth_) {
    iEvent.getByLabel(dtSimLinkTag_,  dtLinksH);
    iEvent.getByLabel(rpcSimLinkTag_, rpcLinksH);
    iEvent.getByLabel(cscSimLinkTag_, cscLinksH);
    iEvent.getByLabel(dtSimHitTag_,   dtSimHitsH);
    iEvent.getByLabel(rpcSimHitTag_,  rpcSimHitsH);
    iEvent.getByLabel(cscSimHitTag_,  cscSimHitsH);
    simHandlesOk = dtLinksH.isValid() && rpcLinksH.isValid() && cscLinksH.isValid();

    if (simHandlesOk) {
      // Build SimTrack.trackId() → gen-muon 1-index map.
      // We consider gen muons with |pdgId|==13 as "gen muons".
      edm::Handle<edm::SimTrackContainer>   simTracksH;
      edm::Handle<reco::GenParticleCollection> genParticlesH;
      iEvent.getByLabel(simTrackTag_,    simTracksH);
      iEvent.getByLabel(genParticleTag_, genParticlesH);

      if (simTracksH.isValid() && genParticlesH.isValid()) {
        // Find the maximum SimTrack trackId to size the vector
        uint32_t maxId = 0;
        for (const auto& st : *simTracksH)
          if (st.trackId() > maxId) maxId = st.trackId();
        simTrackIdToGenIdx.assign(maxId + 1, 0);

        // Build gen-muon list (1-indexed)
        int genMuonIdx = 0;
        for (const auto& gp : *genParticlesH) {
          if (std::abs(gp.pdgId()) != 13)
            continue;
          genMuonIdx++;  // 1-indexed
          // Find all SimTracks that are this gen muon (match via barcode / trackId)
          for (const auto& st : *simTracksH) {
            if (std::abs(st.type()) != 13)
              continue;
            // Match by momentum direction proximity (simple ΔR < 0.05)
            float dEta = st.momentum().eta() - gp.eta();
            float dPhi = st.momentum().phi() - gp.phi();
            // Wrap dPhi
            while (dPhi >  M_PI) dPhi -= 2 * M_PI;
            while (dPhi < -M_PI) dPhi += 2 * M_PI;
            float dR2 = dEta * dEta + dPhi * dPhi;
            if (dR2 < 0.0025f) {  // dR < 0.05
              uint32_t tid = st.trackId();
              if (tid < simTrackIdToGenIdx.size() && simTrackIdToGenIdx[tid] == 0)
                simTrackIdToGenIdx[tid] = genMuonIdx;
            }
          }
        }
      }
    }
  }

  // Flush the per-event cache into OMTFAllInputTree
  for (auto& region : cachedRegions) {
    reg_eventNum   = eventNum;
    reg_iProcessor = region.iProcessor;
    reg_mtfType    = region.mtfType;

    reg_stub_layer   = std::move(region.stub_layer);
    reg_stub_phiHw   = std::move(region.stub_phiHw);
    reg_stub_phiBHw  = std::move(region.stub_phiBHw);
    reg_stub_etaHw   = std::move(region.stub_etaHw);
    reg_stub_r       = std::move(region.stub_r);
    reg_stub_quality = std::move(region.stub_quality);
    reg_stub_type    = std::move(region.stub_type);
    reg_stub_bx      = std::move(region.stub_bx);

    const unsigned int nStubs = reg_stub_layer.size();

    if (doSimTruth_ && simHandlesOk) {
      unsigned int procPhiZero =
          static_cast<unsigned int>(OMTFinputMaker::getProcessorPhiZero(omtfConfig, region.iProcessor));
      assignSimTruth(iEvent, dtLinksH, rpcLinksH, cscLinksH,
                     dtSimHitsH, rpcSimHitsH, cscSimHitsH,
                     simTrackIdToGenIdx, procPhiZero,
                     region.stub_detId, reg_stub_phiHw,
                     reg_stub_trackId, reg_stub_ambiguous,
                     reg_stub_tof, reg_stub_tofSpread, reg_stub_nSimHit);
    } else {
      reg_stub_trackId.assign(nStubs, 0);
      reg_stub_ambiguous.assign(nStubs, 0);
      reg_stub_tof.assign(nStubs, -999.f);
      reg_stub_tofSpread.assign(nStubs, 0.f);
      reg_stub_nSimHit.assign(nStubs, 0);
    }

    allInputTree->Fill();
  }
  cachedRegions.clear();

  // Let the base class write OMTFHitsTree (per-candidate logic unchanged)
  DataROOTDumper2::observeEventEnd(iEvent, finalCandidates);
}

// ---------------------------------------------------------------------------
void DataROOTDumper2AllInput::endJob() {
  allInputTree->Write();
  // Base class writes OMTFHitsTree
  DataROOTDumper2::endJob();
}

