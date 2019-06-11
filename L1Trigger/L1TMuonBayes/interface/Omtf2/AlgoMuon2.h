/*
 * AlgoMuon2.h
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF2_ALGOMUON2_H_
#define INTERFACE_OMTF2_ALGOMUON2_H_

#include <L1Trigger/L1TMuonBayes/interface/AlgoMuonBase.h>
#include "L1Trigger/L1TMuonBayes/interface/MuonStub.h"
#include "L1Trigger/L1TMuonBayes/interface/StubResult.h"

class AlgoMuon2: public AlgoMuonBase {
public:
  AlgoMuon2(const ProcConfigurationBase* config);

  virtual ~AlgoMuon2();

  virtual void addStubResult(float pdfVal, bool valid, int pdfBin, int layer, MuonStubPtr stub);

  int getEtaHw() const override { return etaHw; }

  void setEtaHw() {
    this->etaHw = etaHw;
  }

  bool isValid() const override {
    return valid;
  }

  void setValid(bool valid) {
    this->valid = valid;
  }

  double getPdfSum() const override {
    return pdfSum;
  }

  const bool isKilled() const {
    return killed;
  }

  void kill() {
    killed = true;
    //FIXME maybe also valid = false???
  }


  const StubResult& getStubResult(unsigned int iLayer) const  override {
    return stubResults.at(iLayer);
  }

  const StubResults& getStubResults() const override {
    return stubResults;
  }

  friend std::ostream & operator << (std::ostream &out, const AlgoMuon2& algoMuon2);

  const boost::dynamic_bitset<> getFiredLayerBits() const {
    boost::dynamic_bitset<> firedLayerBitsSum(firedLayerBitsInBx[0].size());
    for(auto& firedLayerBits: firedLayerBitsInBx) {
      firedLayerBitsSum |= firedLayerBits;
    }
    return firedLayerBitsSum;
  }

  int getQuality() const {
    return quality;
  }

  void setQuality(int quality = 0) {
    this->quality = quality;
  }

  int getDisplacementHw() const {
    return displacementHw;
  }

  void setDisplacementHw(int displacementHw = 0) {
    this->displacementHw = displacementHw;
  }

  int getPhiHw() const {
    return phiHw;
  }

  void setPhiHw(int phiHw = 0) {
    this->phiHw = phiHw;
  }

  int getPtHw() const {
    return ptHw;
  }

  void setPtHw(int ptHw = 0) {
    this->ptHw = ptHw;
  }

  int getCharge() const {
    return charge;
  }

  void setCharge(int charge = 0) {
    this->charge = charge;
  }

  double getDisplacedLikelihood() const {
    return displacedLikelihood;
  }

  void setDisplacedLikelihood(double displacedLikelihood = 0) {
    this->displacedLikelihood = displacedLikelihood;
  }

  double getMuonLikelihood() const {
    return muonLikelihood;
  }

  void setMuonLikelihood(double muonLikelihood = 0) {
    this->muonLikelihood = muonLikelihood;
  }

/*  double getSimBeta() const {
    return 0;
  }*/

private:
  int etaHw = 0;

  int phiHw = 0;

  int ptHw = 0;

  int charge = 0;

  int displacementHw = 0;

  bool valid = false;

  double pdfSum = 0;

  double muonLikelihood = 0;
  double displacedLikelihood = 0;

  bool killed = false;

  int quality = 0;

  StubResults stubResults;
};

typedef std::shared_ptr<AlgoMuon2> AlgoMuon2Ptr;
typedef std::vector<AlgoMuon2Ptr> AlgoMuon2s;
#endif /* INTERFACE_OMTF2_ALGOMUON2_H_ */
