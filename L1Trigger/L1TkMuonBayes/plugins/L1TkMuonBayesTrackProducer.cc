#include "L1Trigger/L1TkMuonBayes/plugins/L1TkMuonBayesTrackProducer.h"
#include "L1Trigger/L1TkMuonBayes/interface/PdfModuleWithStats.h"
#include "L1Trigger/L1TkMuonBayes/interface/MuTimingModuleWithStat.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include <iostream>
#include <strstream>
#include <vector>
#include <fstream>
#include <memory>

#include <TFile.h>
#include <TStyle.h>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

L1TkMuonBayesTrackProducer::L1TkMuonBayesTrackProducer(const edm::ParameterSet& cfg)
    : edmParameterSet(cfg), config(std::make_shared<TkMuBayesProcConfig>()) {
  produces<l1t::TkMuonBayesTrackBxCollection>(allTracksProductName);  //all tracks

  produces<l1t::TkMuonBayesTrackBxCollection>(muonTracksProductName);
  //"fast" tracks, i.e. with at least two muon stubs in the same bx as ttRack (i.e. not HSCPs) and passing some cuts

  produces<l1t::TkMuonBayesTrackBxCollection>(hscpTracksProductName);
  //"slow" tracks, i.e. exclusive versus the "fast" tracks and passing some cuts

  MuStubsInputTokens muStubsInputTokens;
  muStubsInputTokens.inputTokenDtPh =
      consumes<L1MuDTChambPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPh"));
  muStubsInputTokens.inputTokenDtTh =
      consumes<L1MuDTChambThContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTTh"));
  muStubsInputTokens.inputTokenCSC =
      consumes<CSCCorrelatedLCTDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcCSC"));
  muStubsInputTokens.inputTokenRPC = consumes<RPCDigiCollection>(edmParameterSet.getParameter<edm::InputTag>("srcRPC"));

  edm::EDGetTokenT<L1Phase2MuDTPhContainer> inputTokenDTPhPhase2;
  if (edmParameterSet.exists("srcDTPhPhase2"))
    inputTokenDTPhPhase2 = consumes<L1Phase2MuDTPhContainer>(edmParameterSet.getParameter<edm::InputTag>("srcDTPhPhase2"));

  edm::InputTag l1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  ttTrackToken = consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >(l1TrackInputTag);

  std::string trackSrc = "L1_TRACKER";
  if (cfg.exists("ttTracksSource"))
    trackSrc = cfg.getParameter<std::string>("ttTracksSource");

  if (trackSrc == "SIM_TRACKS")
    inputTokenSimTracks = mayConsume<edm::SimTrackContainer>(
        edmParameterSet.getParameter<edm::InputTag>("g4SimTrackSrc"));  //TODO is it needed?

  if (trackSrc == "TRACKING_PARTICLES")
    trackingParticleToken = mayConsume<std::vector<TrackingParticle> >(
        edmParameterSet.getParameter<edm::InputTag>("TrackingParticleInputTag"));

  inputMaker = std::make_unique<MuonStubInputMaker>(
      edmParameterSet,
      muStubsInputTokens,
      inputTokenDTPhPhase2,
      config.get(),
      new AngleConverterBase());  //MuonStubInputMaker keeps the AngleConverter

  ttTracksInputMaker = std::make_unique<TTTracksInputMaker>(edmParameterSet);

  //Range of the BXes for which the emulation is performed,
  if (edmParameterSet.exists("bxRangeMin")) {
    bxRangeMin = edmParameterSet.getParameter<int>("bxRangeMin");
  }
  if (edmParameterSet.exists("bxRangeMax")) {
    bxRangeMax = edmParameterSet.getParameter<int>("bxRangeMax");
  }

  if (edmParameterSet.exists("useStubsFromAdditionalBxs")) {
    useStubsFromAdditionalBxs = edmParameterSet.getParameter<int>("useStubsFromAdditionalBxs");
  }

  //muCorrelatorConfig->setBxToProcess(useStubsFromAdditionalBxs + 1); TODO correct, now does not compile due to const

  for (unsigned int ptBin = 0; ptBin < config->getPtHwBins().size(); ++ptBin) {
    edm::LogImportant("l1tOmtfEventPrint") << "ptBin " << setw(2) << ptBin
                                      << " range Hw: " << config->ptBinString(ptBin, 0) << " = "
                                      << config->ptBinString(ptBin, 1) << std::endl;
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TkMuonBayesTrackProducer::~L1TkMuonBayesTrackProducer() {}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TkMuonBayesTrackProducer::beginJob() {
  //m_Reconstruction.beginJob();
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TkMuonBayesTrackProducer::endJob() {
  //m_Reconstruction.endJob();

  IPdfModule* pdfModule = processor->getPdfModule();
  PdfModuleWithStats* pdfModuleWithStats = dynamic_cast<PdfModuleWithStats*>(pdfModule);
  if (pdfModuleWithStats) {
    // using TFileService insteed
    /*gStyle->SetOptStat(111111);
    TFile outfile("muCorrPdfs.root", "RECREATE");
    cout<<__FUNCTION__<<": "<<__LINE__<<" creating file "<<outfile.GetName()<<endl;

    outfile.cd();
    //pdfModuleWithStats->write();*/

    if (edmParameterSet.exists("generatePdfs") && edmParameterSet.getParameter<bool>("generatePdfs")) {
      if (edmParameterSet.exists("outPdfModuleFile")) {
        pdfModuleFile = edmParameterSet.getParameter<std::string>("outPdfModuleFile");
      }
      pdfModuleWithStats->generateCoefficients();
      writePdfs(pdfModule, pdfModuleFile);
    }
  }

  MuTimingModuleWithStat* muTimingModule =
      dynamic_cast<MuTimingModuleWithStat*>(processor->getMuTimingModule());
  if (muTimingModule) {
    if (edmParameterSet.exists("generateTiming") && edmParameterSet.getParameter<bool>("generateTiming")) {
      string muTimingModuleFileName = "muTimingModule.xml";
      if (edmParameterSet.exists(
              "outputTimingFile")) {  //if we read the patterns directly from the xml, we do it only once, at the beginning of the first run, not every run
        muTimingModuleFileName = edmParameterSet.getParameter<std::string>("outputTimingFile");
      }
      muTimingModule->generateCoefficients();
      writeTimingModule(muTimingModule, muTimingModuleFileName);
    }
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TkMuonBayesTrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& eventSetup) {
  if (!processor) {
    std::string pdfModuleType = "PdfModule";
    if (edmParameterSet.exists("pdfModuleType")) {
      pdfModuleType = edmParameterSet.getParameter<std::string>("pdfModuleType");
    }

    PdfModule* pdfModule = nullptr;

    if (pdfModuleType == "PdfModule") {
      pdfModule = new PdfModule(config);
      edm::LogImportant("l1tOmtfEventPrint") << " creating PdfModule for processor" << std::endl;
    } else if (pdfModuleType == "PdfModuleWithStats") {
      edm::LogImportant("l1tOmtfEventPrint") << " creating PdfModuleWithStats for processor" << std::endl;
      pdfModule = new PdfModuleWithStats(config);
    } else {
      throw cms::Exception("L1TkMuonBayesTrackProducer::beginRun: unknown pdfModuleType: " + pdfModuleType);
    }

    if (edmParameterSet.exists("generatePdfs") && edmParameterSet.getParameter<bool>("generatePdfs")) {
      //dont read the pdf if they are going to be generated
    } else if (edmParameterSet.exists("pdfModuleFile")) {
      //if we read the patterns directly from the xml, we do it only once, at the beginning of the first run, not every run
      pdfModuleFile = edmParameterSet.getParameter<edm::FileInPath>("pdfModuleFile").fullPath();
      edm::LogImportant("l1tOmtfEventPrint") << " reading the pdfModule from file " << pdfModuleFile << std::endl;
      readPdfs(pdfModule, pdfModuleFile);
    }

    std::unique_ptr<PdfModule> pdfModuleUniqPtr(pdfModule);

    processor = std::make_unique<TkMuBayesProcessor>(config, std::move(pdfModuleUniqPtr));

    if (edmParameterSet.exists("generateTiming") && edmParameterSet.getParameter<bool>("generateTiming")) {
      std::unique_ptr<MuTimingModule> muTimingModuleUniqPtr =
          std::make_unique<MuTimingModuleWithStat>(config.get());
      processor->setMuTimingModule(muTimingModuleUniqPtr);
    } else if (edmParameterSet.exists("timingModuleFile")) {
      string timingModuleFile = edmParameterSet.getParameter<edm::FileInPath>("timingModuleFile").fullPath();
      edm::LogImportant("l1tOmtfEventPrint")
          << " reading the MuTimingModule from file " << timingModuleFile << std::endl;

      std::unique_ptr<MuTimingModule> muTimingModuleUniqPtr =
          std::make_unique<MuTimingModule>(config.get());
      readTimingModule(muTimingModuleUniqPtr.get(), timingModuleFile);
      processor->setMuTimingModule(muTimingModuleUniqPtr);
    }

    edm::LogImportant("l1tOmtfEventPrint") << " processor constructed" << std::endl;

    //the parameters can be overwritten from the python config
    config->configureFromEdmParameterSet(edmParameterSet);

    inputMaker->initialize(edmParameterSet, eventSetup);
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TkMuonBayesTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup) {
  inputMaker->loadAndFilterDigis(iEvent);

  std::unique_ptr<l1t::TkMuonBayesTrackBxCollection> allTracks(new l1t::TkMuonBayesTrackBxCollection);
  allTracks->setBXRange(bxRangeMin, bxRangeMax);

  std::unique_ptr<l1t::TkMuonBayesTrackBxCollection> muonTracks(new l1t::TkMuonBayesTrackBxCollection);
  muonTracks->setBXRange(bxRangeMin, bxRangeMax);

  std::unique_ptr<l1t::TkMuonBayesTrackBxCollection> hscpTracks(new l1t::TkMuonBayesTrackBxCollection);
  hscpTracks->setBXRange(bxRangeMin, bxRangeMax);

  //std::cout<<"\n"<<__FUNCTION__<<":"<<__LINE__<<" iEvent "<<iEvent.id().event()<<" #####################################################################"<<endl;
  for (int bx = bxRangeMin; bx <= bxRangeMax; bx++) {
    MuonStubsInput muonStubsInput(config.get());
    inputMaker->buildInputForProcessor(
        muonStubsInput.getMuonStubs(), 0, l1t::tftype::bmtf, bx, bx + useStubsFromAdditionalBxs);
    //std::cout<<muonStubsInput<<std::endl;


    LogTrace("l1tOmtfEventPrint") << "\n\nEvent " << iEvent.id().event() << " muonStubsInput bx " << bx << ": \n "
                                  << muonStubsInput << endl;

    auto ttTRacks = ttTracksInputMaker->loadTTTracks(iEvent, bx, edmParameterSet, config.get());
    LogTrace("l1tOmtfEventPrint")<<" ttTRacks.size() " << ttTRacks.size()<< endl;
    for (auto& ttTRack : ttTRacks) {
      LogTrace("l1tOmtfEventPrint") << *ttTRack << endl;
    }
    LogTrace("l1tOmtfEventPrint") << "\n";

    //for(unsigned int iProcessor=0; iProcessor<m_OMTFConfig->nProcessors(); ++iProcessor)
    {
      AlgoTTMuons algoTTMuons = processor->processTracks(muonStubsInput, ttTRacks);
      //fill outgoing collection
      l1t::TkMuonBayesTrackCollection bayesMuCorrTracksInBx = processor->getMuCorrTrackCollection(0, algoTTMuons);
      for (auto& muTrack : bayesMuCorrTracksInBx) {
        allTracks->push_back(bx, muTrack);

        auto firedLayerBits = muTrack.getFiredLayerBits(config->nLayers());
        if (muTrack.hwQual() >= 12 && muTrack.getCandidateType() == l1t::TkMuonBayesTrack::fastTrack &&
            ((firedLayerBits.count() == 2 && muTrack.pdfSum() > 1100) ||
             (firedLayerBits.count() == 3 && muTrack.pdfSum() > 1400) || firedLayerBits.count() >= 4) &&
            ((muTrack.getTtTrackPtr().isNonnull() && muTrack.getTtTrackPtr()->chi2Red() < 200) ||
             muTrack.getTtTrackPtr().isNull())) {
          muonTracks->push_back(bx, muTrack);
        }

        if (muTrack.getCandidateType() == l1t::TkMuonBayesTrack::slowTrack && muTrack.hwQual() >= 13 &&
            ((firedLayerBits.count() == 2 && muTrack.pdfSum() > 1300 && muTrack.getBetaLikelihood() >= 6) ||
             (firedLayerBits.count() == 3 && muTrack.pdfSum() > 1700 && muTrack.getBetaLikelihood() >= 7) ||
             (firedLayerBits.count() == 4 && muTrack.pdfSum() > 2200 && muTrack.getBetaLikelihood() >= 9) ||
             firedLayerBits.count() >= 5) &&
            ((muTrack.getTtTrackPtr().isNonnull() && muTrack.getTtTrackPtr()->chi2Red() < 200) ||
             muTrack.getTtTrackPtr().isNull())
             //todo probably in firmware exactly like that will be not possible, rather cut of chi2 depending on the nStubs
        ) {
          hscpTracks->push_back(bx, muTrack);
        }
      }
    }
  }

  iEvent.put(std::move(allTracks), allTracksProductName);
  iEvent.put(std::move(muonTracks), muonTracksProductName);
  iEvent.put(std::move(hscpTracks), hscpTracksProductName);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void L1TkMuonBayesTrackProducer::readPdfs(IPdfModule* pdfModule, std::string fileName) {
  // open the archive
  std::ifstream ifs(fileName);
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);

  PdfModule* pdfModuleImpl = dynamic_cast<PdfModule*>(pdfModule);
  if (pdfModuleImpl) {
    LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ": " << __LINE__ << " readPdfs from file " << fileName << endl;
    PdfModule& pdfModuleRef = *pdfModuleImpl;
    ia >> BOOST_SERIALIZATION_NVP(pdfModuleRef);

    LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ": " << __LINE__ << " pdfModule->getCoefficients().size() "
                                  << pdfModuleRef.getCoefficients().size() << endl;
  } else {
    throw cms::Exception("L1TkMuonBayesTrackProducer::readPdfs: pdfModule is not of the type PdfModule*");
  }
}

void L1TkMuonBayesTrackProducer::writePdfs(const IPdfModule* pdfModule, std::string fileName) {
  std::ofstream ofs(fileName);

  boost::archive::xml_oarchive xmlOutArch(ofs);
  //boost::archive::text_oarchive txtOutArch(ofs);

  const PdfModule* pdfModuleImpl = dynamic_cast<const PdfModule*>(pdfModule);
  // write class instance to archive
  if (pdfModuleImpl) {
    edm::LogImportant("l1tOmtfEventPrint")
        << __FUNCTION__ << ": " << __LINE__ << " writing pdf to file " << fileName << endl;
    const PdfModule& pdfModuleRef = *pdfModuleImpl;
    xmlOutArch << BOOST_SERIALIZATION_NVP(pdfModuleRef);
    //txtOutArch << (*pdfModuleImpl);
  }
  // archive and stream closed when destructors are called
}

void L1TkMuonBayesTrackProducer::readTimingModule(MuTimingModule* muTimingModule, std::string fileName) {
  // open the archive
  std::ifstream ifs(fileName);
  assert(ifs.good());
  boost::archive::xml_iarchive ia(ifs);

  // write class instance to archive
  //PdfModule* pdfModuleImpl = dynamic_cast<PdfModule*>(pdfModule);
  // write class instance to archive
  if (muTimingModule) {
    LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ": " << __LINE__ << " readTimingModule from file " << fileName
                                  << endl;
    MuTimingModule& muTimingModuleRef = *muTimingModule;
    ia >> BOOST_SERIALIZATION_NVP(muTimingModuleRef);
  } else {
    throw cms::Exception("L1TkMuonBayesTrackProducer::readTimingModule: muTimingModule is 0");
  }
}

void L1TkMuonBayesTrackProducer::writeTimingModule(const MuTimingModule* muTimingModule,
                                                              std::string fileName) {
  std::ofstream ofs(fileName);

  boost::archive::xml_oarchive xmlOutArch(ofs);
  //boost::archive::text_oarchive txtOutArch(ofs);

  //const PdfModule* pdfModuleImpl = dynamic_cast<const PdfModule*>(pdfModule);
  // write class instance to archive
  if (muTimingModule) {
    edm::LogImportant("l1tOmtfEventPrint")
        << __FUNCTION__ << ": " << __LINE__ << " writing MuTimingModule to file " << fileName << endl;
    const MuTimingModule& muTimingModuleRef = *muTimingModule;
    xmlOutArch << BOOST_SERIALIZATION_NVP(muTimingModuleRef);
    //txtOutArch << (*pdfModuleImpl);
  }
  // archive and stream closed when destructors are called
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TkMuonBayesTrackProducer);
