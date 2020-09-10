#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFReconstruction.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/ProcessorBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/XMLConfigReader.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/XMLEventWriter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper2.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EventCapture.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/PatternGenerator.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/PatternOptimizer.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Common/interface/EventBase.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFReconstruction::OMTFReconstruction(const edm::ParameterSet& theConfig, MuStubsInputTokens& muStubsInputTokens) :
  edmParameterSet(theConfig), muStubsInputTokens(muStubsInputTokens), omtfConfig(new OMTFConfiguration()), omtfProc(nullptr), m_OMTFConfigMaker(nullptr) {

  //edmParameterSet.getParameter<std::string>("XMLDumpFileName");
  bxMin = edmParameterSet.exists("bxMin") ? edmParameterSet.getParameter<int>("bxMin") : 0;
  bxMax = edmParameterSet.exists("bxMax") ? edmParameterSet.getParameter<int>("bxMax") : 0;

  inputMaker = std::make_unique<OMTFinputMaker>(theConfig, muStubsInputTokens, omtfConfig.get());
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFReconstruction::~OMTFReconstruction(){

}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::beginJob() {

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::endJob(){
  for(auto& obs : observers) {
    obs->endJob();
  }
}

void OMTFReconstruction::modifyOmtfConfig() {
  if(edmParameterSet.exists("rpcMaxClusterSize") )
    omtfConfig->setRpcMaxClusterSize(edmParameterSet.getParameter<int>("rpcMaxClusterSize"));

  if(edmParameterSet.exists("rpcMaxClusterCnt") )
    omtfConfig->setRpcMaxClusterCnt(edmParameterSet.getParameter<int>("rpcMaxClusterCnt"));

  if(edmParameterSet.exists("rpcDropAllClustersIfMoreThanMax") )
    omtfConfig->setRpcDropAllClustersIfMoreThanMax(edmParameterSet.getParameter<bool>("rpcDropAllClustersIfMoreThanMax"));

  if(edmParameterSet.exists("lctCentralBx")) {
    int lctCentralBx  = edmParameterSet.getParameter<int>("lctCentralBx");
    omtfConfig->setCscLctCentralBx(lctCentralBx);
  }

  if(edmParameterSet.exists("goldenPatternResultFinalizeFunction")) {
    int finalizeFunction = edmParameterSet.getParameter<int>("goldenPatternResultFinalizeFunction");
    omtfConfig->setGoldenPatternResultFinalizeFunction(finalizeFunction);
    edm::LogImportant("OMTFReconstruction") << "GoldenPatternResult::setFinalizeFunction: "<<finalizeFunction << std::endl;
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::beginRun(edm::Run const& run, edm::EventSetup const& eventSetup) {
  const L1TMuonOverlapParams* omtfParams = 0;

  std::string processorType = "OMTFProcessor"; //GoldenPatternWithStat GoldenPattern
  if(edmParameterSet.exists("processorType") ){
    processorType = edmParameterSet.getParameter<std::string>("processorType");
  }

  bool buildPatternsFromXml = (edmParameterSet.exists("patternsXMLFile") || edmParameterSet.exists("patternsXMLFiles"));

  bool firstRun = (omtfProc == 0);

  //if the buildPatternsFromXml == false - we are making the omtfConfig and omtfProc for every run,
  //as the configuration my change between the runs,
  //if buildPatternsFromXml == true - we assume the the entire configuration comes from python,
  //so we do it only for the first run
  if(omtfProc == 0 || buildPatternsFromXml == false) {
    edm::LogImportant("OMTFReconstruction") << "retrieving omtfParams from EventSetup" << std::endl;

    const L1TMuonOverlapParamsRcd& omtfRcd = eventSetup.get<L1TMuonOverlapParamsRcd>();
    edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
    omtfRcd.get(omtfParamsHandle);
    omtfParams = omtfParamsHandle.product();
    if (!omtfParams) {
      edm::LogError("OMTFReconstruction") << "Could not retrieve parameters from Event Setup" << std::endl;
    }
    omtfConfig->configure(omtfParams); //TODO remove, as the configureFromEdmParameterSet does the same

    modifyOmtfConfig(); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //the parameters can be overwritten from the python config
    omtfConfig->configureFromEdmParameterSet(edmParameterSet);

    inputMaker->initialize(edmParameterSet, eventSetup);

    //patterns from the edm::EventSetup are reloaded every beginRun
    if(buildPatternsFromXml == false) {
      edm::LogImportant("OMTFReconstruction") << "getting patterns from EventSetup" << std::endl;
      if(processorType == "OMTFProcessor")
        omtfProc = std::make_unique<OMTFProcessor<GoldenPattern> >(omtfConfig.get(), edmParameterSet, eventSetup, omtfParams);
    }
  }

  if(omtfProc == 0 && buildPatternsFromXml) {//if we read the patterns directly from the xml, we do it only once, at the beginning of the first run, not every run
    std::vector<std::string> patternsXMLFiles;

    if(edmParameterSet.exists("patternsXMLFile")) {
      patternsXMLFiles.push_back(edmParameterSet.getParameter<edm::FileInPath>("patternsXMLFile").fullPath() );
    }
    else if(edmParameterSet.exists("patternsXMLFiles")) {
      for(auto it: edmParameterSet.getParameter<std::vector<edm::ParameterSet> >("patternsXMLFiles")){
        patternsXMLFiles.push_back(it.getParameter<edm::FileInPath>("patternsXMLFile").fullPath());
      }
    }

    for(auto& patternsXMLFile : patternsXMLFiles)
      edm::LogImportant("OMTFReconstruction") << "reading patterns from "<<patternsXMLFile << std::endl;

    XMLConfigReader xmlReader;

    std::string patternType = "GoldenPattern"; //GoldenPatternWithStat GoldenPattern
    if(edmParameterSet.exists("patternType") ){
      patternType = edmParameterSet.getParameter<std::string>("patternType");
    }

    std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    if(patternType == "GoldenPattern") {
      auto gps = xmlReader.readPatterns<GoldenPattern>(*omtfParams, patternsXMLFiles);

      if(processorType == "OMTFProcessor") {
        omtfProc = std::make_unique<OMTFProcessor<GoldenPattern> >(omtfConfig.get(), edmParameterSet, eventSetup, gps);
      }

      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. processorType "<<processorType<<". GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;
    }
    else if(patternType == "GoldenPatternWithStat") {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;

      auto gps = xmlReader.readPatterns<GoldenPatternWithStat>(*omtfParams, patternsXMLFiles);

      if(processorType == "OMTFProcessor") {
        if(edmParameterSet.exists("optimizePatterns") && edmParameterSet.getParameter<bool>("optimizePatterns") ) {
          observers.emplace_back(std::make_unique<PatternOptimizer>(edmParameterSet, omtfConfig.get(), gps));
        }
       
        if(edmParameterSet.exists("generatePatterns") && edmParameterSet.getParameter<bool>("generatePatterns") ) {
          observers.emplace_back(std::make_unique<PatternGenerator>(edmParameterSet, omtfConfig.get(), gps));
          edm::LogImportant("OMTFReconstruction") << "generatePatterns: true " << std::endl;
        }
        
        omtfProc = std::make_unique<OMTFProcessor<GoldenPatternWithStat> >(omtfConfig.get(), edmParameterSet, eventSetup, gps) ;
      }
    }
    else if(patternType == "GoldenPatternWithThresh") {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
      auto gps = xmlReader.readPatterns<GoldenPatternWithThresh>(*omtfParams, patternsXMLFiles);
      
      omtfProc = std::make_unique<OMTFProcessor<GoldenPatternWithThresh> >(omtfConfig.get(), edmParameterSet, eventSetup, gps) ;
      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;
    }
    else {
      throw cms::Exception("OMTFReconstruction::beginRun: unknown GoldenPattern type: " + patternType);
    }
  }

  omtfProc->printInfo();

  if(firstRun) {
    addObservers();
  }

}

void OMTFReconstruction::addObservers() {
  if(edmParameterSet.exists("dumpResultToXML")) {
    if(edmParameterSet.getParameter<bool>("dumpResultToXML"))
      observers.emplace_back(std::make_unique<XMLEventWriter>(omtfConfig.get(), edmParameterSet.getParameter<std::string>("XMLDumpFileName") ));
  }

  if(edmParameterSet.exists("eventCaptureDebug"))
    if(edmParameterSet.getParameter<bool>("eventCaptureDebug") )
      observers.emplace_back(std::make_unique<EventCapture>(edmParameterSet, omtfConfig.get() ) );

  if( edmParameterSet.exists("dumpResultToROOT"))
    if(edmParameterSet.getParameter<bool>("dumpResultToROOT"))
      observers.emplace_back(std::make_unique<DataROOTDumper>(edmParameterSet, omtfConfig.get() ));

  auto omtfProcGoldenPat = dynamic_cast<OMTFProcessor<GoldenPattern>* >(omtfProc.get());
  if(omtfProcGoldenPat) {
    if(edmParameterSet.exists("dumpHitsToROOT") && edmParameterSet.getParameter<bool>("dumpHitsToROOT")) {
      std::string rootFileName = edmParameterSet.getParameter<std::string>("dumpHitsFileName");
      observers.emplace_back(std::make_unique<DataROOTDumper2>(edmParameterSet, omtfConfig.get(), omtfProcGoldenPat->getPatterns(), rootFileName) );
    }
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::unique_ptr<l1t::RegionalMuonCandBxCollection> OMTFReconstruction::reconstruct(const edm::Event& iEvent, const edm::EventSetup& evSetup) {
  LogTrace("l1tMuBayesEventPrint")<<"\n"<<__FUNCTION__<<":"<<__LINE__<<" iEvent "<<iEvent.id().event()<<endl;
  inputMaker->loadAndFilterDigis(iEvent);
  
  for(auto& obs : observers) {
    obs->observeEventBegin(iEvent);
  }

  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates = std::make_unique<l1t::RegionalMuonCandBxCollection>();
  candidates->setBXRange(bxMin, bxMax);


  ///The order is important: first put omtf_pos candidates, then omtf_neg.
  for(int bx = bxMin; bx<= bxMax; bx++) {
  
    for(unsigned int iProcessor=0; iProcessor<omtfConfig->nProcessors(); ++iProcessor) {
      std::vector<l1t::RegionalMuonCand> candMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_pos, bx, inputMaker.get(), observers);

      //fill outgoing collection
      for (auto & candMuon :  candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    for(unsigned int iProcessor=0; iProcessor<omtfConfig->nProcessors(); ++iProcessor) {
      std::vector<l1t::RegionalMuonCand> candMuons = omtfProc->run(iProcessor, l1t::tftype::omtf_neg, bx, inputMaker.get(), observers);

      //fill outgoing collection
      for (auto & candMuon :  candMuons) {
        candidates->push_back(bx, candMuon);
      }
    }

    //edm::LogInfo("OMTFReconstruction") <<"OMTF:  Number of candidates in BX="<<bx<<": "<<candidates->size(bx) << std::endl;;
  }
  
  for(auto& obs : observers) {
    obs->observeEventEnd(iEvent, candidates);
  }

  return candidates;
}


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

