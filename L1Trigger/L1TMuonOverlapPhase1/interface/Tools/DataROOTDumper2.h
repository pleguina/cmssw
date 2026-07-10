/*
 * DataROOTDumper2.h
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#ifndef L1T_OmtfP1_TOOLS_DATAROOTDUMPER2_H_
#define L1T_OmtfP1_TOOLS_DATAROOTDUMPER2_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EmulationObserverBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/CandidateSimMuonMatcher.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "TMap.h"
#include "TArrayI.h"
#include "TFile.h"
#include "TH2.h"

#include <functional>

class TTree;

struct OmtfEvent {
public:
  unsigned int eventNum = 0;

  //muonPt = 0 means that no muon was matched to the candidate
  short muonEvent = -1;

  //Explicit "is this a gen-matched entry?" flag. All muon* / vertex* / deltaEta / deltaPhi fields above
  //are optional gen-truth-match fields: when hasGenMatch == false they hold their "no match" sentinel
  //values (muonEvent=-1, muonPt=0, etc.) and the row still corresponds to a real OMTF candidate (or, in
  //samples with no truth source configured at all, to every OMTF candidate in the event). This makes
  //fake/ghost studies possible on samples where, by construction, nothing is gen-matched (e.g. pure
  //background/noise productions) since such candidates are no longer silently dropped.
  bool hasGenMatch = false;

  float muonPt = 0, muonEta = 0, muonPhi = 0, muonPropEta = 0, muonPropPhi = 0;
  char muonCharge = 0;
  float muonDxy = 0;
  float muonRho = 0;
  //pdgId in principle should be int
  short parentPdgId = 0;
  float vertexEta = 0;
  float vertexPhi = 0;

  float omtfPt = 0, omtfEta = 0, omtfPhi = 0, omtfUPt = 0;
  char omtfCharge = 0;
  char omtfProcessor = 0;
  short omtfScore = 0;
  short omtfRefHitPhi = 0;

  short omtfHwEta = 0;

  char omtfQuality = 0;
  char omtfRefLayer = 0;
  char omtfRefHitNum = 0;

  unsigned int omtfFiredLayers = 0;

  bool killed = false;

  float deltaPhi = 0, deltaEta = 0;

  //float omtfPtCont = 0;

  struct Hit {
    union {
      unsigned long rawData = 0;

      struct {
        char layer;
        char quality;
        char etaHw;
        char valid;
        short deltaR;  // stub_r - refStub_r  [cm]
        short phiDist; // stub_phi - refStub_phi  [HW units]; phiBHw for bending layers
      };
    };

    ~Hit() {}
  };

  std::vector<unsigned long> hits;

  // Extended per-hit branches — parallel to 'hits' (same index i → same stub)
  // These carry the fields that do not fit in the 64-bit Hit union.
  std::vector<short> hits_phiBHw;  // DT bending angle in HW units; 0 for non-DT layers
  std::vector<short> hits_phiHw;   // absolute stub phi in HW units (before ref-subtraction)
  std::vector<short> hits_r;       // absolute radial distance [cm]
  std::vector<signed char> hits_type;  // MuonStub::Type enum cast to int8
  std::vector<signed char> hits_bx;    // bunch crossing offset from BX=0
};

class DataROOTDumper2 : public EmulationObserverBase {
public:
  DataROOTDumper2(const edm::ParameterSet& edmCfg,
                  const OMTFConfiguration* omtfConfig,
                  CandidateSimMuonMatcher* candidateSimMuonMatcher);

  ~DataROOTDumper2() override;

  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>&,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const FinalMuons& finalMuons) override;

  void observeEventEnd(const edm::Event& iEvent,
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override;

  void endJob() override;

private:
  void initializeTTree();

  CandidateSimMuonMatcher* candidateSimMuonMatcher = nullptr;

  TTree* rootTree = nullptr;

  OmtfEvent omtfEvent;

  unsigned int evntCnt = 0;

  TH1I* ptGenPos = nullptr;
  TH1I* ptGenNeg = nullptr;

  //std::vector<TH2*> hitVsPt;

  bool dumpKilledOmtfCands = false;
};

#endif /* L1T_OmtfP1_TOOLS_DATAROOTDUMPER2_H_ */
