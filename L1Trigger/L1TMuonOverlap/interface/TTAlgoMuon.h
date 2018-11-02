/*
 * TTAlgoMuon.h
 *
 *  Created on: Oct 19, 2018
 *      Author: kbunkow
 */

#ifndef OMTF_TTALGOMUON_H_
#define OMTF_TTALGOMUON_H_

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"

class TrackingTriggerTrack {
public:
/*  TrackingTriggerTrack(double phi, double eta,  double pt, int charge) : phi(phi), eta(eta), pt(pt), charge(charge) {
    //todo convert to hw scale
  };

  TrackingTriggerTrack(int phiHw, int etaHw,  int ptHw, int charge) : charge(charge), phiHw(phiHw), etaHw(etaHw), ptHw(ptHw) {
    //todo convert to physics scale
  };*/

  TrackingTriggerTrack(const SimTrack& simMuon);

  TrackingTriggerTrack(const TTTrack< Ref_Phase2TrackerDigi_>& ttTRack, int l1Tk_nPar);

  int getCharge() const {
    return charge;
  }

  double getEta() const {
    return eta;
  }

  double getPhi() const {
    return phi;
  }

  double getPt() const {
    return pt;
  }

  int getEtaHw() const {
    return etaHw;
  }

  int getPhiHw() const {
    return phiHw;
  }

  int getPtHw() const {
    return ptHw;
  }

  void setEtaHw(int etaHw = 0) {
    this->etaHw = etaHw;
  }

  void setPhiHw(int phiHw = 0) {
    this->phiHw = phiHw;
  }

  void setPtHw(int ptHw = 0) {
    this->ptHw = ptHw;
  }

private:
  double phi = 0;
  double eta = 0;
  double pt = 0;
  int charge = 0;


  ///in integer hardware scales
  int phiHw = 0;
  int etaHw = 0;
  int ptHw = 0;

  int indx = -1;
};

typedef std::vector<TrackingTriggerTrack> TTTracks;

class TTAlgoMuon: public AlgoMuon {
public:
  //move the gpResults content to the this->gpResults, gpResults is empty after that
  TTAlgoMuon(const TrackingTriggerTrack& ttTrack, const GoldenPatternResult& gpResult, GoldenPatternBase* goldenPatern,
      std::vector<std::shared_ptr<GoldenPatternResult> >& gpResults, unsigned int refHitNum, int bx = 0):
    AlgoMuon(gpResult, goldenPatern,  refHitNum, bx),
    ttTrack(ttTrack) {
    this->gpResults.swap(gpResults);
  }

  const TrackingTriggerTrack& getTtTrack() const {
    return ttTrack;
  }

  ///
  void setGpResults(std::vector<std::shared_ptr<GoldenPatternResult> >& gpResults) {
    this->gpResults.swap(gpResults);
    gpResults.clear(); //in case there was something before in the this->gpResults
  }

  std::vector<std::shared_ptr<GoldenPatternResult> >& getGpResults() {
    return gpResults;
  }

  //TODO add getters for phi eta pt at vertex?

private:
  TrackingTriggerTrack ttTrack; //maybe rather should be a pointer?

  ///results for all the reference layers processed for this ttTrack
  ///we cannot use the GoldenPatternResuls stored by the GoldenPatternBase, because in one event the same goldePattern
  ///may be hit many times by different ttTracks
  std::vector<std::shared_ptr<GoldenPatternResult> > gpResults;
};

typedef std::vector<std::shared_ptr<TTAlgoMuon> > TTMuons;


#endif /* OMTF_TTALGOMUON_H_ */
