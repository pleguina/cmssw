/*
 * AlgoMuonBase.h
 *
 *  Created on: Mar 1, 2019
 *      Author: Karol Bunkowski kbunkow@cern.ch
 */

#ifndef INTERFACE_ALGOMUONBASE_H_
#define INTERFACE_ALGOMUONBASE_H_

#include "L1Trigger/L1TMuonBayes/interface/MuonStub.h"
#include "L1Trigger/L1TMuonBayes/interface/StubResult.h"

class AlgoMuonBase {
public:
  AlgoMuonBase();
  virtual ~AlgoMuonBase();

  virtual bool isValid() const = 0;

  //virtual void setValid(bool valid) = 0;

  virtual unsigned int getFiredLayerCnt() const = 0;

  virtual double getPdfSum() const = 0;

 /* virtual const bool isKilled() const = 0;

  virtual void kill()  = 0;*/

  //virtual bool isLayerFired(unsigned int iLayer) const  = 0;

  virtual const StubResult& getStubResult(unsigned int iLayer) const  = 0;

  virtual const StubResults& getStubResults() const  = 0;

  //virtual const boost::dynamic_bitset<>& getFiredLayerBits() const  = 0;

};

#endif /* INTERFACE_ALGOMUONBASE_H_ */
