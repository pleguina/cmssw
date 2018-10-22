#ifndef OMTFReconstruction_H
#define OMTFReconstruction_H

#include "xercesc/util/XercesDefs.hpp"

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

#include "L1Trigger/L1TMuonOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/GhostBuster.h"
#include "L1Trigger/L1TMuonOverlap/interface/IProcessorEmulator.h"
#include "L1Trigger/L1TMuonOverlap/interface/IOMTFEmulationObserver.h"

class L1TMuonOverlapParams;
class OMTFConfiguration;
class OMTFConfigMaker;
class XMLConfigWriter;

namespace XERCES_CPP_NAMESPACE{
  class DOMElement;
  class DOMDocument;
  class DOMImplementation;
}

class OMTFReconstruction {
  public:
    //OMTFReconstruction();

    OMTFReconstruction(const edm::ParameterSet&);

    ~OMTFReconstruction();

    void beginJob();

    void endJob();

    void beginRun(edm::Run const& run, edm::EventSetup const& iSetup);  

    std::unique_ptr<l1t::RegionalMuonCandBxCollection> reconstruct(const edm::Event&, const edm::EventSetup&);

  private:

    edm::ParameterSet m_Config;

/*    void loadAndFilterDigis(const edm::Event&);

    void getProcessorCandidates(unsigned int iProcessor, l1t::tftype mtfType, int bx,
            l1t::RegionalMuonCandBxCollection & myCandidates);*/
  

    bool dumpResultToXML, dumpDetailedResultToXML;
    int bxMin, bxMax;

  ///OMTF objects
    OMTFConfiguration   *m_OMTFConfig;

    //OMTFProcessor<GoldenPattern>  *m_OMTF;
    unique_ptr<IProcessorEmulator> m_OMTF;
  ///
    //xercesc::DOMElement *aTopElement;
    OMTFConfigMaker     *m_OMTFConfigMaker;
    //XMLConfigWriter     *m_Writer;

    std::vector<std::unique_ptr<IOMTFEmulationObserver> > observers;
    unsigned int theEvent;
};

#endif
