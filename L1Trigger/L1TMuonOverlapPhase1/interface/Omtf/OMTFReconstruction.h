#ifndef OMTFReconstruction_H
#define OMTFReconstruction_H

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GhostBuster.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IProcessorEmulator.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFProcessor.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

class OMTFConfiguration;
class OMTFConfigMaker;

class OMTFReconstruction {
  public:

    OMTFReconstruction(const edm::ParameterSet&, MuStubsInputTokens& muStubsInputTokens);

    virtual ~OMTFReconstruction();

    void beginJob();

    void endJob();

    void beginRun(edm::Run const& run, edm::EventSetup const& iSetup);  

    std::unique_ptr<l1t::RegionalMuonCandBxCollection> reconstruct(const edm::Event&, const edm::EventSetup&);

    void virtual modifyOmtfConfig();

    //takes the ownership of the inputMaker
    void setInputMaker(OMTFinputMaker* inputMaker) {
      this->inputMaker.reset(inputMaker);
    }

    void virtual addObservers();
  protected:

    edm::ParameterSet edmParameterSet;

    MuStubsInputTokens& muStubsInputTokens;


    int bxMin, bxMax;

  ///OMTF objects
    unique_ptr<OMTFConfiguration> omtfConfig;

    unique_ptr<OMTFinputMaker> inputMaker;

    unique_ptr<IProcessorEmulator> omtfProc;

    OMTFConfigMaker* m_OMTFConfigMaker;

    std::vector<std::unique_ptr<IOMTFEmulationObserver> > observers;
};

#endif
