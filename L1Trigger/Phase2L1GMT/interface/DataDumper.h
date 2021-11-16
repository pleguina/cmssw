/*
 * DataDumper.h
 *
 *  Created on: Nov 5, 2021
 *      Author: kbunkow
 */

#ifndef INTERFACE_DATADUMPER_H_
#define INTERFACE_DATADUMPER_H_

#include "L1Trigger/Phase2L1GMT/interface/PreTrackMatchedMuon.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"

#include "TTree.h"
#include "TFile.h"

namespace Phase2L1GMT {

const unsigned int tfLayersCnt = 5;
struct TrackMatchedMuonRecord {
  //from trackingParticle
  float tpPt = 0, tpEta = 0, tpPhi = 0;
  short tpCharge = 0;
  char type = 0;

  //0 - no match, 1 - very loose, 2 - loose , 3
  char matching = 0;

  //from ttTrack
  uint tttCharge = 0;
  uint tttPt = 0;
  int tttEta = 0;
  int tttPhi = 0;
  int tttZ0 = 0;
  int tttD0 = 0;

  //from GMT
  uint beta = 0;
  bool isGlobal = 0;
  uint quality = 0;

  std::vector<unsigned char> deltaCoords1;
  std::vector<unsigned char> deltaCoords2;

  std::vector<char> stubTiming;
  //std::vector<unsigned char> stubTiming2;

  TrackMatchedMuonRecord(): deltaCoords1(tfLayersCnt), deltaCoords2(tfLayersCnt), stubTiming(tfLayersCnt) {}
};


class PreTrackMatchedMuonProcessor {
public:
  PreTrackMatchedMuonProcessor() {};
  virtual ~PreTrackMatchedMuonProcessor() {};

  virtual void process(PreTrackMatchedMuon& preTrackMatchedMuon) = 0;
};

class DataDumper: public PreTrackMatchedMuonProcessor {
public:
  DataDumper(const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& ttTrackMCTruthToken,
             const edm::EDGetTokenT< std::vector< TrackingParticle > >& trackingParticleToken, std::string rootFileName);

  virtual ~DataDumper();

  void getHandles(const edm::Event& event);

  void process(PreTrackMatchedMuon& preTrackMatchedMuon) override;

private:
  void initializeTTree(std::string rootFileName);
  void saveTTree();

  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken;
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > mcTruthTTTrackHandle;

  edm::EDGetTokenT< std::vector< TrackingParticle > > trackingParticleToken;
  edm::Handle< std::vector< TrackingParticle > > trackingParticleHandle;

  std::vector<edm::Ptr< TrackingParticle > > muonTrackingParticles;
  bool muonTrackingParticlesFilled = false;

  TFile* rootFile = nullptr;
  TTree* rootTree = nullptr;

  TrackMatchedMuonRecord record;

  unsigned int evntCnt = 0;

  //TH1I* ptGenPos = nullptr;
  //TH1I* ptGenNeg = nullptr;
};

} /* namespace Phase2L1GMT */

#endif /* INTERFACE_DATADUMPER_H_ */
