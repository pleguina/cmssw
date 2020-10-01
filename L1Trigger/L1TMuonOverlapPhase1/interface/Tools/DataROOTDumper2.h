/*
 * DataROOTDumper2.h
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_TOOLS_DATAROOTDUMPER2_H_
#define INTERFACE_TOOLS_DATAROOTDUMPER2_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/GpResultsToPt.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/PatternOptimizerBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "TMap.h"
#include "TArrayI.h"

#include <functional>

class TTree;

struct OmtfEvent {

public:

  double muonPt = 0, muonEta = 0, muonPhi = 0;
  int muonCharge = 0;

  int omtfCharge = 0, omtfProcessor = 0, omtfScore = 0;
  double omtfPt = 0, omtfEta = 0, omtfPhi = 0;
  unsigned int omtfQuality = 0, omtfRefLayer = 0;
  unsigned int omtfFiredLayers = 0;

  float omtfPtCont = 0;

  struct Hit {
    union {
      unsigned long rawData;

      struct {
        char layer = 0;
        char quality = 0;
        char z = 0;
        char valid = 0;
        short eta = 0;
        short phiDist = 0;
      };
    };

  };

  std::vector<unsigned long> hits;

  boost::multi_array<float, 2> omtfGpResultsPdfSum; //[iREfHit][iGp]
  boost::multi_array<unsigned int, 2> omtfGpResultsFiredLayers; //[iREfHit][iGp]

  OmtfEvent(unsigned int nRefHits, unsigned int nGoldenPatterns):
    omtfGpResultsPdfSum(boost::extents[nRefHits][nGoldenPatterns]),
    omtfGpResultsFiredLayers(boost::extents[nRefHits][nGoldenPatterns]) {}
};

class DataROOTDumper2: public PatternOptimizerBase {
public:
  DataROOTDumper2(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
         const std::vector<std::shared_ptr<GoldenPattern> >& gps, std::string rootFileName);

  virtual ~DataROOTDumper2();

  virtual void observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates);

  virtual void endJob();

private:
  std::vector<std::shared_ptr<GoldenPattern> > gps;

  void initializeTTree(std::string rootFileName);
  void saveTTree();

  TFile* rootFile = nullptr;
  TTree* rootTree = nullptr;

  OmtfEvent event;

  unsigned int evntCnt = 0;

  TH1I* ptGenPos = nullptr;
  TH1I* ptGenNeg = nullptr;

  bool dumpGpResults = false;

  GpResultsToPt* gpResultsToPt = nullptr; //TODO move to OmtfProcessor
};

#endif /* INTERFACE_TOOLS_DATAROOTDUMPER2_H_ */
