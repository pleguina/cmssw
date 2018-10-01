#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/plugins/L1TMuonOverlapTrackProducer.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFReconstruction.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"
#include "L1Trigger/L1TMuonOverlap/interface/OmtfName.h"
#include "L1Trigger/L1TMuonOverlap/interface/GhostBusterPreferRefDt.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorterWithThreshold.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigReader.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternParametrised.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlap/interface/PatternOptimizer.h"

#include "L1Trigger/L1TMuonOverlap/interface/XMLEventWriter.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include <boost/timer/timer.hpp>

/*OMTFReconstruction::OMTFReconstruction() :
  m_OMTFConfig(0), m_OMTF(0), m_OMTFConfigMaker(0){}*/
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFReconstruction::OMTFReconstruction(const edm::ParameterSet& theConfig) :
  m_Config(theConfig), m_OMTFConfig(0), m_OMTF(0), m_OMTFConfigMaker(0) {

  dumpResultToXML = m_Config.getParameter<bool>("dumpResultToXML");
  dumpDetailedResultToXML = m_Config.getParameter<bool>("dumpDetailedResultToXML");
  m_Config.getParameter<std::string>("XMLDumpFileName");  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFReconstruction::~OMTFReconstruction(){
  
  delete m_OMTFConfig;
  delete m_OMTF;  

  //if (m_Writer) delete m_Writer;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::beginJob() {
    //std::cout<<__FUNCTION__<<":"<<__LINE__<<"test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    m_OMTFConfig = new OMTFConfiguration();
    //m_OMTF = new OMTFProcessor<GoldenPattern>();
    //m_OMTF = new OMTFProcessor<GoldenPatternParametrised>();

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::endJob(){
/*  if(dumpResultToXML){
    std::string fName = m_Config.getParameter<std::string>("XMLDumpFileName");
    m_Writer->finaliseXMLDocument(fName);
  } */

  for(auto& obs : observers) {
    obs->endJob();
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
  const L1TMuonOverlapParams* omtfParams = 0;
  if(m_OMTF == 0 || m_Config.exists("patternsXMLFile") == false) {
    edm::LogImportant("OMTFReconstruction") << "retrieving parameters from Event Setup" << std::endl;

    const L1TMuonOverlapParamsRcd& omtfRcd = iSetup.get<L1TMuonOverlapParamsRcd>();
    edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
    omtfRcd.get(omtfParamsHandle);
    omtfParams = omtfParamsHandle.product();
    if (!omtfParams) {
      edm::LogError("OMTFReconstruction") << "Could not retrieve parameters from Event Setup" << std::endl;
    }
    m_OMTFConfig->configure(omtfParams);

    //patterns from the L1TMuonOverlapParamsESProducer, are reloaded every run begin
    if(m_OMTF == 0 && m_Config.exists("patternsXMLFile") == false) {
      edm::LogInfo("OMTFReconstruction") << "getting patterns from L1TMuonOverlapParamsESProducer" << std::endl;
      //m_OMTFConfig->initPatternPtRange();
      OMTFProcessor<GoldenPattern>* proc =  new OMTFProcessor<GoldenPattern>(m_OMTFConfig);
      m_OMTF = proc;

      m_OMTF->configure(m_OMTFConfig, omtfParams);
      m_OMTFConfig->setPatternPtRange(proc->getPatternPtRange() );
    }
  }
  if(m_OMTF == 0 && m_Config.exists("patternsXMLFile") ) {//if we read the patterns directly from the xml, we do it only once, at the beginning of the first run, not every run
    std::string patternsXMLFile = m_Config.getParameter<edm::FileInPath>("patternsXMLFile").fullPath();
    edm::LogImportant("OMTFReconstruction") << "reading patterns from "<<patternsXMLFile << std::endl;
    XMLConfigReader xmlReader;
    xmlReader.setPatternsFile(patternsXMLFile);

    std::string patternType = "GoldenPattern"; //GoldenPatternParametrised GoldenPatternWithStat GoldenPattern
    if(m_Config.exists("patternType") ){
      patternType = m_Config.getParameter<std::string>("patternType");
    }
    std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    if(patternType == "GoldenPattern") {
      OMTFProcessor<GoldenPattern>* proc =  new OMTFProcessor<GoldenPattern>(m_OMTFConfig);
      m_OMTF = proc;

      auto const& gps = xmlReader.readPatterns<GoldenPattern>(*omtfParams);
      proc->setGPs(gps);
      m_OMTFConfig->setPatternPtRange(proc->getPatternPtRange() );
      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;
    }
    else if(patternType == "GoldenPatternWithStat") {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
      OMTFProcessor<GoldenPatternWithStat>* proc =  new OMTFProcessor<GoldenPatternWithStat>(m_OMTFConfig);
      m_OMTF = proc;

      auto gps = xmlReader.readPatterns<GoldenPatternWithStat>(*omtfParams);
      proc->setGPs(gps);
      m_OMTFConfig->setPatternPtRange(proc->getPatternPtRange() );
      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;

      std::unique_ptr<IOMTFEmulationObserver> obs(new PatternOptimizer(m_Config, m_OMTFConfig, gps));
      observers.emplace_back(std::move(obs));

      for(auto gp : gps) {
        gp->init();
      }

      for(auto& gp : gps) {
        edm::LogImportant("OMTFReconstruction")<<gp->key()<<" "
            <<m_OMTFConfig->getPatternPtRange(gp->key().theNumber).ptFrom
            <<" - "<<m_OMTFConfig->getPatternPtRange(gp->key().theNumber).ptTo<<" GeV"<<std::endl;
      }

      if(m_Config.exists("sorterType") ) {//TODO add it also for the patternType == "GoldenPattern" - if needed
        string sorterType = m_Config.getParameter<std::string>("sorterType");
        edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. sorterType: "<<sorterType<< std::endl;
        if(sorterType == "sorterWithThreshold") {
          GoldenPatternResult::setFinalizeFunction(2);

          OMTFSorterWithThreshold<GoldenPatternWithStat>::Mode mode = OMTFSorterWithThreshold<GoldenPatternWithStat>::bestGPByMaxGpProbability1;
          string modeStr = m_Config.getParameter<std::string>("sorterWithThresholdMode");
          if(modeStr == "bestGPByThresholdOnProbability2")
            mode = OMTFSorterWithThreshold<GoldenPatternWithStat>::bestGPByThresholdOnProbability2;
          else if(modeStr == "bestGPByMaxGpProbability1")
            mode = OMTFSorterWithThreshold<GoldenPatternWithStat>::bestGPByMaxGpProbability1;

          proc->setSorter(new OMTFSorterWithThreshold<GoldenPatternWithStat>(m_OMTFConfig, mode));
        }
      }
    }
    else if(patternType == "GoldenPatternWithThresh") {
      std::cout<<__FUNCTION__<<":"<<__LINE__<<std::endl;
      OMTFProcessor<GoldenPatternWithThresh>* proc =  new OMTFProcessor<GoldenPatternWithThresh>(m_OMTFConfig);
      m_OMTF = proc;

      auto gps = xmlReader.readPatterns<GoldenPatternWithThresh>(*omtfParams);
      proc->setGPs(gps);
      m_OMTFConfig->setPatternPtRange(proc->getPatternPtRange() );
      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;

      //std::unique_ptr<IOMTFEmulationObserver> obs(new PatternOptimizer(m_Config, m_OMTFConfig, gps));
      //observers.emplace_back(std::move(obs));

      for(auto& gp : gps) {
        edm::LogImportant("OMTFReconstruction")<<gp->key()<<" "
            <<m_OMTFConfig->getPatternPtRange(gp->key().theNumber).ptFrom
            <<" - "<<m_OMTFConfig->getPatternPtRange(gp->key().theNumber).ptTo<<" GeV"<<" threshold "<<gp->getThreshold(0)<<std::endl;
      }

      if(m_Config.exists("sorterType") ) {//TODO add it also for the patternType == "GoldenPattern" - if needed
        string sorterType = m_Config.getParameter<std::string>("sorterType");
        edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. sorterType: "<<sorterType<< std::endl;
        if(sorterType == "sorterWithThreshold") {
          GoldenPatternResult::setFinalizeFunction(2);

          OMTFSorterWithThreshold<GoldenPatternWithThresh>::Mode mode = OMTFSorterWithThreshold<GoldenPatternWithThresh>::bestGPByMaxGpProbability1;
          string modeStr = m_Config.getParameter<std::string>("sorterWithThresholdMode");
          if(modeStr == "bestGPByThresholdOnProbability2")
            mode = OMTFSorterWithThreshold<GoldenPatternWithThresh>::bestGPByThresholdOnProbability2;
          else if(modeStr == "bestGPByMaxGpProbability1")
            mode = OMTFSorterWithThreshold<GoldenPatternWithThresh>::bestGPByMaxGpProbability1;

          proc->setSorter(new OMTFSorterWithThreshold<GoldenPatternWithThresh>(m_OMTFConfig, mode));
        }
      }
    }
    else if(patternType == "GoldenPatternParametrised") {
      OMTFProcessor<GoldenPatternParametrised>* proc = new OMTFProcessor<GoldenPatternParametrised>(m_OMTFConfig);
      m_OMTF = proc;

      auto const& gps = xmlReader.readPatterns<GoldenPattern>(*omtfParams);
      proc->resetConfiguration();
      for(auto& gp :  gps) {
        if(gp.get() != 0) {
          gp->setConfig(m_OMTFConfig);
          edm::LogImportant("OMTFReconstruction") <<gp->key()<< std::endl;
          GoldenPatternParametrised* newGp = new GoldenPatternParametrised(gp.get());
          proc->addGP(newGp);
        }
      }

      m_OMTFConfig->setPatternPtRange(proc->getPatternPtRange() );
      edm::LogImportant("OMTFReconstruction") << "OMTFProcessor constructed. GoldenPattern type: "<<patternType<<" size: "<<gps.size() << std::endl;
    }
    else {
      throw cms::Exception("OMTFReconstruction::beginRun: unknown GoldenPattern type: " + patternType);
    }

  }

  m_InputMaker.initialize(iSetup, m_OMTFConfig);

  if(dumpResultToXML){
    std::unique_ptr<IOMTFEmulationObserver> obs(new XMLEventWriter(m_OMTFConfig, m_Config.getParameter<std::string>("XMLDumpFileName")));
    observers.emplace_back(std::move(obs));
  }

  if(m_Config.exists("ghostBusterType") ) {
    if(m_Config.getParameter<std::string>("ghostBusterType") == "GhostBusterPreferRefDt")
      m_OMTF->setGhostBuster(new GhostBusterPreferRefDt(m_OMTFConfig));
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::unique_ptr<l1t::RegionalMuonCandBxCollection> OMTFReconstruction::reconstruct(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

  loadAndFilterDigis(iEvent);

  //if(dumpResultToXML) aTopElement = m_Writer->writeEventHeader(iEvent.id().event());

  for(auto& obs : observers) {
    obs->observeEventBegin(iEvent);
  }

  // NOTE: assuming all is for bx 0
  int bx = 0;
  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates(new l1t::RegionalMuonCandBxCollection);

  ///The order is important: first put omtf_pos candidates, then omtf_neg.
  for(unsigned int iProcessor=0; iProcessor<m_OMTFConfig->nProcessors(); ++iProcessor)
    getProcessorCandidates(iProcessor, l1t::tftype::omtf_pos, bx, *candidates);

  for(unsigned int iProcessor=0; iProcessor<m_OMTFConfig->nProcessors(); ++iProcessor)
    getProcessorCandidates(iProcessor, l1t::tftype::omtf_neg, bx, *candidates);
    
  for(auto& obs : observers) {
    obs->observeEventEnd(iEvent);
  }

  return candidates;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::loadAndFilterDigis(const edm::Event& iEvent){

  // Filter digis by dropping digis from selected (by cfg.py) subsystems
  if(!m_Config.getParameter<bool>("dropDTPrimitives")){
    iEvent.getByLabel(m_Config.getParameter<edm::InputTag>("srcDTPh"),dtPhDigis);
    iEvent.getByLabel(m_Config.getParameter<edm::InputTag>("srcDTTh"),dtThDigis);
  }
  if(!m_Config.getParameter<bool>("dropRPCPrimitives")) iEvent.getByLabel(m_Config.getParameter<edm::InputTag>("srcRPC"),rpcDigis);
  if(!m_Config.getParameter<bool>("dropCSCPrimitives")) iEvent.getByLabel(m_Config.getParameter<edm::InputTag>("srcCSC"),cscDigis);
  
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFReconstruction::getProcessorCandidates(unsigned int iProcessor, l1t::tftype mtfType, int bx,
               l1t::RegionalMuonCandBxCollection & omtfCandidates){

  //boost::timer::auto_cpu_timer t("%ws wall, %us user in getProcessorCandidates\n");
  m_InputMaker.setFlag(0);
  OMTFinput input = m_InputMaker.buildInputForProcessor(dtPhDigis.product(),
                dtThDigis.product(),
                cscDigis.product(),
                rpcDigis.product(),
                iProcessor, mtfType);
  int flag = m_InputMaker.getFlag();
  //cout<<"buildInputForProce "; t.report();
  m_OMTF->processInput(iProcessor, mtfType, input);
  //cout<<"processInput       "; t.report();
  std::vector<AlgoMuon> algoCandidates =  m_OMTF->sortResults(iProcessor, mtfType);
  //cout<<"sortResults        "; t.report();
  // perform GB 
  std::vector<AlgoMuon> gbCandidates =  m_OMTF->ghostBust(algoCandidates);
  //cout<<"ghostBust          "; t.report();
  // fill RegionalMuonCand colleciton
  std::vector<l1t::RegionalMuonCand> candMuons = m_OMTF->getFinalcandidates(iProcessor, mtfType, gbCandidates);
  //cout<<"getFinalcandidates "; t.report();
  //fill outgoing collection
  for (auto & candMuon :  candMuons) {
     candMuon.setHwQual( candMuon.hwQual() | flag);         //FIXME temporary debug fix
     omtfCandidates.push_back(bx, candMuon);
  }
  //dump to XML
  //writeResultToXML(iProcessor, mtfType,  input, algoCandidates, candMuons);
  for(auto& obs : observers) {
    obs->observeProcesorEmulation(iProcessor, mtfType,  input, algoCandidates, gbCandidates, candMuons);
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
