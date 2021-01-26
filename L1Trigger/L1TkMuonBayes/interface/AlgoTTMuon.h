/*
 * AlgoTTMuon.h
 *
 *  Created on: Feb 1, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef L1TkMuonBayes_ALGOTTMUON_H_
#define L1TkMuonBayes_ALGOTTMUON_H_

#include "L1Trigger/L1TkMuonBayes/interface/TkMuBayesProcConfig.h"
#include "L1Trigger/L1TkMuonBayes/interface/TrackingTriggerTrack.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/AlgoMuonBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/StubResult.h"

#include <vector>
#include <memory>
#include "boost/dynamic_bitset.hpp"

class AlgoTTMuon : public AlgoMuonBase {
public:
  AlgoTTMuon(const TrackingTriggerTrackPtr& ttTrack, TkMuBayesProcConfigPtr& config)
      : AlgoMuonBase(config.get()), ttTrack(ttTrack), stubResults(config->nLayers()){};

  AlgoTTMuon(const TrackingTriggerTrackPtr& ttTrack, TkMuBayesProcConfigPtr& config, const MuonStubPtr& refStub)
      : AlgoMuonBase(config.get()), ttTrack(ttTrack), stubResults(config->nLayers()), refStub(refStub){};

  ~AlgoTTMuon() override{};

  virtual void addStubResult(float pdfVal, bool valid, int pdfBin, int layer, MuonStubPtr stub);

  virtual void invalidateStubResult(int layer);

  int getEtaHw() const override { return ttTrack->getEtaHw(); }

  bool isValid() const override { return valid; }

  void setValid(bool valid) { this->valid = valid; }

  double getPdfSum() const override { return pdfSum; }

  const bool isKilled() const { return killed; }

  void kill() {
    killed = true;
    //FIXME maybe also valid = false???
  }

  const TrackingTriggerTrackPtr& getTTTrack() const { return ttTrack; }

  const StubResult& getStubResult(unsigned int iLayer) const override { return stubResults.at(iLayer); }

  const StubResults& getStubResults() const override { return stubResults; }

  friend std::ostream& operator<<(std::ostream& out, const AlgoTTMuon& algoTTMuon);

  int getQuality() const { return quality; }

  void setQuality(int quality = 0) { this->quality = quality; }

  double getSimBeta() const override { return ttTrack->getSimBeta(); }

private:
  bool valid = false;

  double pdfSum = 0;

  bool killed = false;

  int quality = 0;

  ///Number of fired layers - excluding bending layers
  //unsigned int firedLayerCnt = 0;

  //ttTrack, stubResults and refStub should be needed in the emulation (debugging etc), but not in the firmware

  TrackingTriggerTrackPtr ttTrack;

  StubResults stubResults;

  MuonStubPtr refStub;
};

typedef std::shared_ptr<AlgoTTMuon> AlgoTTMuonPtr;
typedef std::vector<AlgoTTMuonPtr> AlgoTTMuons;

#endif /* L1TkMuonBayes_ALGOTTMUON_H_ */
