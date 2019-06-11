/*
 * Omtf2inputMaker.h
 *
 *  Created on: May 15, 2019
 *      Author: kbunkow
 */

#ifndef INTERFACE_OMTF2_OMTF2INPUTMAKER_H_
#define INTERFACE_OMTF2_OMTF2INPUTMAKER_H_

#include "L1Trigger/L1TMuonBayes/interface/Omtf/OMTFinputMaker.h"

class Omtf2InputMaker: public OMTFinputMaker {
public:
  Omtf2InputMaker();

  virtual ~Omtf2InputMaker();

  virtual void addStub(MuonStubPtrs2D& muonStubsInLayers, unsigned int iLayer, unsigned int iInput, MuonStub& stub) override;
};

#endif /* INTERFACE_OMTF2_OMTF2INPUTMAKER_H_ */
