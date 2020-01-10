/*
 * DataROOTDumper2.h
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTFPATTERNGENERATION_DATAROOTDUMPER2_H_
#define INTERFACE_OMTFPATTERNGENERATION_DATAROOTDUMPER2_H_



#include "L1Trigger/L1TMuonBayes/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonBayes/interface/OmtfPatternGeneration/PatternOptimizerBase.h"
#include <functional>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "TMap.h"
#include "TArrayI.h"

class TTree;

struct OmtfEvent {

public:

  double muonPt, muonEta, muonPhi;
  int muonCharge;

  int omtfCharge, omtfProcessor, omtfScore;
  double omtfPt, omtfEta, omtfPhi;
  unsigned int omtfQuality, omtfRefLayer, omtfHitsWord;


  struct Hit {
    union {
      unsigned long rawData;

      struct {
        char layer = 0;
        char quality = 0;
        char z = 0;
        short eta = 0;
        short phiDist = 0;
      };
    };

  };

  std::vector<unsigned long> hits;

};

class DataROOTDumper2: public PatternOptimizerBase {
public:
  DataROOTDumper2(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig,
         std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps, std::string rootFileName);

  virtual ~DataROOTDumper2();

  virtual void observeEventEnd(const edm::Event& iEvent, std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates);

  virtual void endJob();

private:

  void initializeTTree(std::string rootFileName);
  void saveTTree();

  TFile* rootFile = nullptr;
  TTree* rootTree = nullptr;

  OmtfEvent event;

  unsigned int evntCnt = 0;

  TH1I* ptGenPos = nullptr;
  TH1I* ptGenNeg = nullptr;
};

#endif /* INTERFACE_OMTFPATTERNGENERATION_DATAROOTDUMPER2_H_ */
