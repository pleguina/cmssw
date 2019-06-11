#include <iostream>
#include <strstream>
#include <vector>
#include <fstream>
#include <memory>

#include <TFile.h>
#include <TStyle.h>

#include "Math/VectorUtil.h"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFinput.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFProcessor.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/XMLConfigWriter.h>
#include <L1Trigger/L1TMuonBayes/plugins/L1TMuonBayesOmtf2TrackProducer.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

L1TMuonBayesOmtf2TrackProducer::L1TMuonBayesOmtf2TrackProducer(const edm::ParameterSet& cfg)
  : edmParameterSet(cfg),
    omtfConfig(std::make_unique<OMTFConfiguration>()),
    inputMaker(std::make_unique<Omtf2InputMaker>())
{
  produces<l1t::RegionalMuonCandBxCollection >("OMTF2");

  muStubsInputTokens.inputTokenDTPh = consumes<L1MuDTChambPhContainer>(cfg.getParameter<edm::InputTag>("srcDTPh"));
  muStubsInputTokens.inputTokenDTTh = consumes<L1MuDTChambThContainer>(cfg.getParameter<edm::InputTag>("srcDTTh"));
  muStubsInputTokens.inputTokenCSC = consumes<CSCCorrelatedLCTDigiCollection>(cfg.getParameter<edm::InputTag>("srcCSC"));
  muStubsInputTokens.inputTokenRPC = consumes<RPCDigiCollection>(cfg.getParameter<edm::InputTag>("srcRPC"));

  genParticleToken = consumes<reco::GenParticleCollection>(edmParameterSet.getParameter<edm::InputTag>("genParticle"));
  simTrackToken = consumes<edm::SimTrackContainer>(edmParameterSet.getParameter<edm::InputTag>("g4SimTrackSrc")); //TODO remove

  if(cfg.exists("etaCutFrom") && cfg.exists("etaCutTo")) {
    etaCutFrom = cfg.getParameter< double >("etaCutFrom");
    etaCutTo = cfg.getParameter< double >("etaCutTo");
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TMuonBayesOmtf2TrackProducer::~L1TMuonBayesOmtf2TrackProducer(){
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtf2TrackProducer::beginJob(){



}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtf2TrackProducer::endJob(){
  if(ptModuleLut2DGen) {
    ptModuleLut2DGen->endJob();
  }

  processor->getPtModule()->saveAsRootHists();


  std::string sampleTypeStr = "muons";
  if(sampleType == PtModuleLut2DGen::SampleType::displacedMuons)
    sampleTypeStr = "displacedMuon";

  std::string modeStr = "GenMean";
  if(mode == PtModuleLut2DGen::Mode::calcualteStdDevs)
    modeStr = "GenStdDevs";

  std::string fileName = "ptModuleLut2D_" + modeStr + "_" +  sampleTypeStr + ".xml";

  if(fileName != "")
    writePtModule(processor->getPtModule().get(), fileName);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtf2TrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& evSetup){
  //inputMaker = std::make_unique<Omtf2inputMaker>();

  if(processor) //at the moment do only once per job
    return;

  const L1TMuonOverlapParams* omtfParams = 0;

  const L1TMuonOverlapParamsRcd& omtfRcd = evSetup.get<L1TMuonOverlapParamsRcd>();
  edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
  omtfRcd.get(omtfParamsHandle);
  omtfParams = omtfParamsHandle.product();
  if (!omtfParams) {
    edm::LogError("OMTFReconstruction") << "Could not retrieve parameters from Event Setup" << std::endl;
  }
  omtfConfig->configure(omtfParams);

  inputMaker->initialize(edmParameterSet, evSetup, this->omtfConfig.get(), muStubsInputTokens);


  std::string ptModuleType = "PtModuleLut2D"; //GoldenPatternParametrised GoldenPatternWithStat GoldenPattern
  if(edmParameterSet.exists("ptModuleType") ) {
    ptModuleType = edmParameterSet.getParameter<std::string>("ptModuleType");
  }

  if(ptModuleType == "PtModuleLut2DGen") {
    ptModuleLut2DGen = std::make_shared<PtModuleLut2DGen>(omtfConfig.get());
    processor = std::make_unique<OmtfProcessorLut2d>(omtfConfig.get(), ptModuleLut2DGen);
  }
  else
    processor = std::make_unique<OmtfProcessorLut2d>(omtfConfig.get());

  edm::LogImportant("l1tMuBayesEventPrint")<<" processor constructed with "<<ptModuleType<<std::endl;

  std::string ptModuleFile = "";

  if(edmParameterSet.exists("ptModuleFile") ) {
    ptModuleFile = edmParameterSet.getParameter<edm::FileInPath>("ptModuleFile").fullPath();
    edm::LogImportant("l1tMuBayesEventPrint")<<" reading the ptModule from file "<<ptModuleFile<<std::endl;
    readPtModule(processor->getPtModule().get(), ptModuleFile);
  }
  else {
    edm::LogImportant("l1tMuBayesEventPrint")<<" no ptModule file given. Calling  processor->getPtModule()->init()"<<std::endl;
    processor->getPtModule()->init();
  }

  edm::LogImportant("l1tMuBayesEventPrint")<<" processor->getPtModule()->getLut2ds().size() = "<<processor->getPtModule()->getLut2ds().size()<<std::endl;
  if(processor->getPtModule()->getLut2ds().size())
    edm::LogImportant("l1tMuBayesEventPrint")<<" LUT output count = "<<processor->getPtModule()->getLut2ds().at(0)->getLutValues().size()<<std::endl;

  if(ptModuleLut2DGen) {
    std::string genModeStr = edmParameterSet.getParameter<std::string>("genMode");
    std::string sampleTypeStr = edmParameterSet.getParameter<std::string>("sampleType");

    if(genModeStr == "calcualteMeans")
      mode = PtModuleLut2DGen::Mode::calcualteMeans;
    else if(genModeStr == "calcualteStdDevs")
      mode = PtModuleLut2DGen::Mode::calcualteStdDevs;

    if(sampleTypeStr == "muons")
      sampleType = PtModuleLut2DGen::SampleType::muons;
    else if(sampleTypeStr == "displacedMuon")
      sampleType = PtModuleLut2DGen::SampleType::displacedMuons;

    ptModuleLut2DGen->prepare(mode, sampleType);
  }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonBayesOmtf2TrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup) {
  std::ostringstream str;
  inputMaker->loadAndFilterDigis(iEvent);

  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates(new l1t::RegionalMuonCandBxCollection);
  candidates->setBXRange(bxRangeMin, bxRangeMax);

  //std::cout<<"\n"<<__FUNCTION__<<":"<<__LINE__<<" iEvent "<<iEvent.id().event()<<" #####################################################################"<<endl;

  LogTrace("l1tMuBayesEventPrint")<<"\n\nEvent "<<iEvent.id().event()<<endl;
  if(ptModuleLut2DGen) {
    /*std::vector<const SimTrack*> smiMuons = findSimMuon(iEvent);
    if(smiMuons.size() > 0) {
      ptModuleLut2DGen->setSimMuon(smiMuons[0]); //todo what to do if more then 1 simMuon?
    }
    else {
      ptModuleLut2DGen->setSimMuon(nullptr);
    } */

    std::vector<const reco::GenParticle*> genMuons = findGenMuon(iEvent);
    if(genMuons.size() > 0) {
      ptModuleLut2DGen->setGenMuon(genMuons[0]); //todo what to do if more then 1 simMuon?
    }
    else {
      ptModuleLut2DGen->setGenMuon(nullptr);
    }
  }


  for(int bx = bxRangeMin; bx <= bxRangeMax; bx++) {
    //std::cout<<muonStubsInput<<std::endl;
    for(int side = 0; side <= 1; side++)  {
      for(unsigned int iProcessor=0; iProcessor<omtfConfig->nProcessors(); ++iProcessor) {
        MuonStubsInput muonStubsInput(omtfConfig.get());
        l1t::tftype procType = side == 0 ? l1t::tftype::omtf_neg : l1t::tftype::omtf_pos;
        inputMaker->buildInputForProcessor(muonStubsInput.getMuonStubs(), iProcessor, procType, bx, bx + useStubsFromAdditionalBxs);

        bool areStubs = false;
        for(auto& layer : muonStubsInput.getMuonStubs()) {
          if(layer.size() )
            areStubs = true;
        } //TODO crate function in MuonStubsInput

        if(areStubs) {
          LogTrace("l1tMuBayesEventPrint")<<"iProcessor "<<iProcessor<<" muonStubsInput in bx "<<bx<<": \n "<<muonStubsInput<<endl;
          LogTrace("l1tMuBayesEventPrint")<<"\n";
        }

        AlgoMuon2s algoMuons = processor->processStubs(muonStubsInput);

        std::vector<l1t::RegionalMuonCand> candMuons = processor->getFinalcandidates(iProcessor, procType, algoMuons);

        if(areStubs)
          LogTrace("l1tMuBayesEventPrint")<<"\n";

        //fill outgoing collection
        for (auto & candMuon :  candMuons) {
          candidates->push_back(bx, candMuon);
        }
      }
    }
  }

  //TODO
//  std::unique_ptr<l1t::RegionalMuonCandBxCollection > candidates = m_Reconstruction.reconstruct(iEvent, evSetup);
  iEvent.put(std::move(candidates), "OMTF2");
}

void L1TMuonBayesOmtf2TrackProducer::readPtModule(PtModuleLut2D* ptModule, std::string fileName) {
  // open the archive
  std::ifstream ifs(fileName);
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);

  // write class instance to archive
  //PdfModule* pdfModuleImpl = dynamic_cast<PdfModule*>(pdfModule);
  // write class instance to archive
  if(ptModule) {
    LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<" readPtModule from file "<<fileName<<endl;
    PtModuleLut2D& ptModuleRef = *ptModule;
    ia >> BOOST_SERIALIZATION_NVP(ptModuleRef);
  }
  else {
    throw cms::Exception("L1TMuonBayesMuCorrelatorTrackProducer::readPtModule: ptModule is 0");
  }

  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<" ptModule red successfully, ptModule->getLut2ds().size() "<<ptModule->getLut2ds().size()<<endl;
  boost::dynamic_bitset<> interpolate(ptModule->getLut2ds().at(0)->getLutConfig().outValCnt, 0xffffff); //TODO change if needed
  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<std::endl;
  ptModule->setInterpolate(interpolate);
  LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<std::endl;
}

void L1TMuonBayesOmtf2TrackProducer::writePtModule(const PtModuleLut2D* ptModule, std::string fileName) {
  std::ofstream ofs(fileName);

  boost::archive::xml_oarchive xmlOutArch(ofs);
  //boost::archive::text_oarchive txtOutArch(ofs);

  //const PdfModule* pdfModuleImpl = dynamic_cast<const PdfModule*>(pdfModule);
  // write class instance to archive
  if(ptModule) {
    edm::LogImportant("l1tMuBayesEventPrint")<<__FUNCTION__<<": "<<__LINE__<<" writing ptModule to file "<<fileName<<endl;
    const PtModuleLut2D& ptModuleRef = *ptModule;
    xmlOutArch << BOOST_SERIALIZATION_NVP(ptModuleRef);
    //txtOutArch << (*pdfModuleImpl);
  }
  // archive and stream closed when destructors are called
}

//const std::vector<const SimTrack*> L1TMuonBayesOmtf2TrackProducer::findSimMuon(const edm::Event &event) {
//  std::vector<const SimTrack*> muons;
//
//  if(edmParameterSet.exists("g4SimTrackSrc") == false)
//    return muons;
//
//  edm::Handle<edm::SimTrackContainer> simTks;
//  event.getByToken(simTrackToken, simTks);
//
//  for (std::vector<SimTrack>::const_iterator it=simTks->begin(); it< simTks->end(); it++) {
//    const SimTrack* simMu = nullptr;
//
//    const SimTrack& aTrack = *it;
//    if ( !(aTrack.type() == 13 || aTrack.type() == -13) )
//      continue;
//
//    simMu = &aTrack;
//
//    bool wasCloseMuon = false;
//    for(auto& previous : muons) {
//      if(previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(), previous->momentum()) < 0.07) { //TODO what sense it has?
//        previous = nullptr;
//        wasCloseMuon = true;
//        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" removing previous muon with pt "<<previous->momentum().pt()
//            <<" beacuse there is anoter close muon with pt "<<simMu->momentum().pt()<<std::endl;
//      /*if (previous.momentum().pt() > result->momentum().pt())
//        wasCloseMuon = true;*/
//      }
//    }
//
//    //aTrack.
//
//    int genPartIndx = simMu->genpartIndex();
//    if(genPartIndx >= 0) {
//      /*LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticles->size() "<<genParticles->size()
//          <<" simTks->size() "<<simTks->size()<<" genPartIndx "<<genPartIndx<<std::endl;
//      const reco::GenParticle& genPart = (*genParticles)[genPartIndx];
//
//      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle:"
//          <<" px "<<genPart.px()<<" py "<<genPart.py()<<" pt "<<genPart.pt()
//          <<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<endl;*/
//
//      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" simTrack   : pdgId "<<simMu->type()
//          <<" px "<<simMu->momentum().Px()<<" py "<<simMu->momentum().Py()<<" pt "<<simMu->momentum().pt()<<endl;
//    }
//    if(!wasCloseMuon)
//      muons.push_back(simMu);
//  }
//
///*  for(size_t i = 0; i < genParticles->size() ; ++ i) {
//      const reco::GenParticle& genPart = (*genParticles)[i];
//      //int id = p.pdgId();
//      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
//          <<" px "<<genPart.px()<<" py "<<genPart.py()<<" pz "<<genPart.pz()<<" pt "<<genPart.pt()
//          <<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<endl;
//
//      LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
//          <<" px "<<genPart.px()/genPart.p()<<" py "<<genPart.py()/genPart.p()<<" pz "<<genPart.pz()/genPart.p()<<" pt "<<genPart.pt()<<endl;
//          //<<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<endl;
//  }*/
//  return muons;
//}

const std::vector<const reco::GenParticle*> L1TMuonBayesOmtf2TrackProducer::findGenMuon(const edm::Event &event) {
  std::vector<const reco::GenParticle*> muons;

  if(edmParameterSet.exists("genParticle") == false)
    return muons;

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(genParticleToken, genParticles);


  for(size_t i = 0; i < genParticles->size() ; ++ i) {
    const reco::GenParticle& genPart = (*genParticles)[i];

    if( abs(genPart.pdgId() )== 13 ) {
      if( abs (genPart.momentum().eta() ) >= etaCutFrom && abs (genPart.momentum().eta() ) <= etaCutTo) {
        genPart.momentum().eta();

        muons.push_back( &genPart);
        //int id = p.pdgId();
        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
                  <<" px "<<genPart.px()<<" py "<<genPart.py()<<" pz "<<genPart.pz()<<" pt "<<genPart.pt()
                  <<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<" eta "<<genPart.momentum().eta()<<endl;

        LogTrace("l1tMuBayesEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
                  <<" px "<<genPart.px()/genPart.p()<<" py "<<genPart.py()/genPart.p()<<" pz "<<genPart.pz()/genPart.p()<<" pt "<<genPart.pt()<<endl;
        //<<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<endl;
      }
    }
  }

  return muons;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonBayesOmtf2TrackProducer);
