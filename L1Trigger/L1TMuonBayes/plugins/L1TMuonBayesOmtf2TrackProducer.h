#ifndef OMTFProducer_H
#define OMTFProducer_H

#include <L1Trigger/L1TMuonBayes/interface/Omtf2/Omtf2InputMaker.h>
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/OmtfProcessorLut2d.h"
#include "L1Trigger/L1TMuonBayes/interface/Omtf2/PtModuleLut2DGen.h"

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

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class L1TMuonOverlapParams;
class OMTFConfiguration;

class L1TMuonBayesOmtf2TrackProducer : public edm::EDProducer {
 public:
  L1TMuonBayesOmtf2TrackProducer(const edm::ParameterSet&);

  ~L1TMuonBayesOmtf2TrackProducer() override;

  void beginJob() override;

  void endJob() override;

  void beginRun(edm::Run const& run, edm::EventSetup const& iSetup) override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  void readPtModule(PtModuleLut2D* muTimingModule, std::string fileName);

  void writePtModule(const PtModuleLut2D* muTimingModule, std::string fileName);

  //const std::vector<const SimTrack*> findSimMuon(const edm::Event &event);

  const std::vector<const reco::GenParticle*> findGenMuon(const edm::Event &event);

 private:
  edm::ParameterSet edmParameterSet;
  
  MuStubsInputTokens muStubsInputTokens;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken;

  //bool dumpResultToXML, dumpDetailedResultToXML;

  std::unique_ptr<OMTFConfiguration> omtfConfig;

  std::unique_ptr<Omtf2InputMaker> inputMaker;
  std::unique_ptr<OmtfProcessorLut2d> processor;

  std::shared_ptr<PtModuleLut2DGen> ptModuleLut2DGen;

  //Range of the BXes for which the emulation is performed,
  int bxRangeMin = 0, bxRangeMax = 0;

  //if 1 then the emulator takes the input data from one more BX, which allows to reconstruct the HSCPs
  int useStubsFromAdditionalBxs = 0;

  std::string ptModuleType;

  PtModuleLut2DGen::Mode mode = PtModuleLut2DGen::Mode::calcualteMeans;
  PtModuleLut2DGen::SampleType sampleType = PtModuleLut2DGen::SampleType::muons;

  double etaCutFrom = 0.8;
  double etaCutTo = 1.24;
};

#endif
