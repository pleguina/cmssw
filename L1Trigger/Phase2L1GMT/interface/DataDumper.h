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

struct TrackMatchedMuonRecord {
  //from trackingParticle
  unsigned int eventNum = 0;
  short tpEvent = -1;
  float tpPt = 0, tpEta = 0, tpPhi = 0, tpBeta = 0;
  short tpCharge = 0;
  int tpType = 0;

  //0 - no match, 1 - very loose, 2 - loose , 3
  char matching = 0;

  //from ttTrack
  short tttCurvature = 0;
  char tttCharge = 0;
  unsigned short tttPt = 0;
  short tttEta = 0;
  short tttPhi = 0;
  short tttZ0 = 0;
  short tttD0 = 0;

  //from TTTrack_TrackWord
  char tttChi2rphi = 0; //4 bits
  char tttChi2rz = 0; //4 bits
  char tttBendChi2 = 0; //3 bits
  char tttQualityMVA = 0; //3 bits
  char tttOtherMVA = 0; //6 bits

  //from GMT
  uint gmtBeta = 0;
  bool isGlobal = 0;
  uint quality = 0;

  unsigned short hitsValid = 0;

  std::vector<unsigned char> deltaCoords1;
  std::vector<unsigned char> deltaCoords2;

  std::vector<unsigned char> deltaEtas1;
  std::vector<unsigned char> deltaEtas2;
  
  std::vector<char> stubTiming;
  std::vector<char> stubType;


  //std::vector<unsigned char> stubTiming2;

  TrackMatchedMuonRecord():
    deltaCoords1(TF_LAYER_CNT), deltaCoords2(TF_LAYER_CNT),
    deltaEtas1(TF_LAYER_CNT), deltaEtas2(TF_LAYER_CNT),
    stubTiming(TF_LAYER_CNT), stubType(TF_LAYER_CNT) {}
};


class PreTrackMatchedMuonProcessor {
public:
  PreTrackMatchedMuonProcessor() {};
  virtual ~PreTrackMatchedMuonProcessor() {};

  virtual void eventBegin() = 0;

  virtual void process(PreTrackMatchedMuon& preTrackMatchedMuon) = 0;
};

class DataDumper: public PreTrackMatchedMuonProcessor {
public:
  DataDumper(const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >& ttTrackMCTruthToken,
             const edm::EDGetTokenT< std::vector< TrackingParticle > >& trackingParticleToken, std::string rootFileName);

  virtual ~DataDumper();

  void getHandles(const edm::Event& event);

  void eventBegin() override {
    evntCnt++;
  }

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
