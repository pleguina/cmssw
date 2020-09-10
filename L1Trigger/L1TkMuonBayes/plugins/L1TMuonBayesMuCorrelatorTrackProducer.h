#ifndef L1TMuonOverlapTTMergerTrackProducer_H
#define L1TMuonOverlapTTMergerTrackProducer_H

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "L1Trigger/L1TkMuonBayes/interface/MuCorrelatorConfig.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStubMakerBase.h"
#include "L1Trigger/L1TkMuonBayes/interface/MuCorrelatorInputMaker.h"
#include "L1Trigger/L1TkMuonBayes/interface/MuCorrelatorProcessor.h"
#include "L1Trigger/L1TkMuonBayes/interface/TTTracksInputMaker.h"

#include <memory>
#include <string>
#include <vector>

class L1TMuonBayesMuCorrelatorTrackProducer : public edm::EDProducer {
 public:
  L1TMuonBayesMuCorrelatorTrackProducer(const edm::ParameterSet&);

  ~L1TMuonBayesMuCorrelatorTrackProducer() override;

  void beginJob() override;

  void endJob() override;

  void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  static constexpr char allTracksProductName[] = "AllTracks"; //all tracks produced by the muon correlator, without additional cuts
  static constexpr char muonTracksProductName[] = "MuonTracks"; //"fast" tracks, i.e. with at least two muon stubs in the same bx as ttRack (=> not HSCPs) and with some cuts reducing rate
  static constexpr char hscpTracksProductName[] = "HscpTracks"; //"slow" tracks, i.e. exclusive versus the "fast" tracks and passing some cuts

 private:


  void readPdfs(IPdfModule* pdfModule, std::string fileName);
  void writePdfs(const IPdfModule* pdfModule, std::string fileName);

  void readTimingModule(MuTimingModule* muTimingModule, std::string fileName);
  void writeTimingModule(const MuTimingModule* muTimingModule, std::string fileName);

  edm::ParameterSet edmParameterSet;
  
/*  edm::EDGetTokenT<L1MuDTChambPhContainer> inputTokenDTPh;
  edm::EDGetTokenT<L1MuDTChambThContainer> inputTokenDTTh;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCSC;
  edm::EDGetTokenT<RPCDigiCollection> inputTokenRPC;*/

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken;

  edm::EDGetTokenT<edm::SimTrackContainer> inputTokenSimTracks; //TODO remove

  edm::EDGetTokenT< std::vector< TrackingParticle > > trackingParticleToken;

  bool dumpResultToXML = false;

  MuCorrelatorConfigPtr muCorrelatorConfig;

  std::unique_ptr<MuCorrelatorInputMaker> inputMaker;
  std::unique_ptr<MuCorrelatorProcessor> muCorrelatorProcessor;
  std::unique_ptr<TTTracksInputMaker> ttTracksInputMaker;

  std::string pdfModuleFile = "pdfModule.xml";

  //Range of the BXes for which the emulation is performed,
  int bxRangeMin = 0, bxRangeMax = 0;

  //if 1 then the emulator takes the input data from one more BX, which allows to reconstruct the HSCPs
  int useStubsFromAdditionalBxs = 0;
};

#endif

