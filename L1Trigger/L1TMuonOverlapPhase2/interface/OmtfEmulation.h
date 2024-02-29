/*
 * OmtfEmulation.h
 *
 *  Created on: May 20, 2020
 *      Author: kbunkow
 */

#ifndef L1Trigger_L1TMuonOverlapPhase2_OmtfEmulation_h
#define L1Trigger_L1TMuonOverlapPhase2_OmtfEmulation_h

#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFReconstruction.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfPhase2AngleConverter.h"
#include "L1Trigger/L1TMuonOverlapPhase2/interface/InputMakerPhase2.h"

class OmtfEmulation : public OMTFReconstruction {
public:
  OmtfEmulation(const edm::ParameterSet& edmParameterSet,
                MuStubsInputTokens& muStubsInputTokens,
                MuStubsPhase2InputTokens& muStubsPhase2InputTokens);

  void beginJob();

  ~OmtfEmulation() override = default;

  void addObservers(const MuonGeometryTokens& muonGeometryTokens,
                    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>& magneticFieldEsToken,
                    const edm::ESGetToken<Propagator, TrackingComponentsRecord>& propagatorEsToken) override;

private:
  MuStubsPhase2InputTokens& muStubsPhase2InputTokens;
  unique_ptr<PtAssignmentBase> ptAssignment;
};

#endif /* L1Trigger_L1TMuonOverlapPhase2_OmtfEmulation_h */
