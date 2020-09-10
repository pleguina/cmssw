/*
 * OmtfEmulation.h
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#ifndef L1TMUONOVERLAPPHASE1_OMTFEMULATION_H_
#define L1TMUONOVERLAPPHASE1_OMTFEMULATION_H_

#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFReconstruction.h"

class OmtfEmulation: public OMTFReconstruction {
public:
	OmtfEmulation(const edm::ParameterSet& edmParameterSet, MuStubsInputTokens& muStubsInputTokens, edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2);

	virtual ~OmtfEmulation();

	void addObservers();
};

#endif /* L1TMUONOVERLAPPHASE1_OMTFEMULATION_H_ */
