/*
 * EmulationObserverBase.cc
 *
 *  Created on: Aug 18, 2021
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/EmulationObserverBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Math/VectorUtil.h"

EmulationObserverBase::EmulationObserverBase(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig)
    : edmCfg(edmCfg), omtfConfig(omtfConfig), simMuon(nullptr) {}

EmulationObserverBase::~EmulationObserverBase() {
  // TODO Auto-generated destructor stub
}

void EmulationObserverBase::observeProcesorEmulation(unsigned int iProcessor,
                                                     l1t::tftype mtfType,
                                                     const std::shared_ptr<OMTFinput>& input,
                                                     const AlgoMuons& algoCandidates,
                                                     const AlgoMuons& gbCandidates,
                                                     const std::vector<l1t::RegionalMuonCand>& candMuons) {
  unsigned int procIndx = omtfConfig->getProcIndx(iProcessor, mtfType);

  /*
  double ptSim = simMuon->momentum().pt();
  int chargeSim = (abs(simMuon->type()) == 13) ? simMuon->type()/-13 : 0;
  int patNum = omtfConfig->getPatternNum(ptSim, chargeSim);
  GoldenPatternWithStat* exptCandGp = goldenPatterns.at(patNum).get(); // expected pattern
*/

  //bool found = false;

  unsigned int i = 0;
  for (auto& gbCandidate : gbCandidates) {
    //int iRefHit = gbCandidate.getRefHitNumber();
    if (gbCandidate->getGoldenPatern() != nullptr &&
        gbCandidate->getGpResult().getFiredLayerCnt() > omtfCand->getGpResult().getFiredLayerCnt()) {
      //edm::LogVerbatim("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" gbCandidate "<<gbCandidate<<" "<<std::endl;
      omtfCand = gbCandidate;
      //omtfResult = gbCandidate.getGoldenPatern()->getResults()[procIndx][iRefHit]; //TODO be carrefful, because in principle the results sored by the goldenPattern can be altered in one event. In phae I omtf this should not happened, but in OMTFProcessorTTMerger - yes
      //exptResult = exptCandGp->getResults()[procIndx][iRefHit];
      candProcIndx = procIndx;

      regionalMuonCand = candMuons.at(i);
      //should be good, as the regionalMuonCand is created for every  gbCandidate in OMTFProcessor<GoldenPatternType>::getFinalcandidates
      //found = true;

      //this->algoCandidates = algoCandidates; //TODO uncomment if needed
    }
    i++;
  }

  //////////////////////debug printout/////////////////////////////
  /*if(found) {
    GoldenPatternWithStat* omtfCandGp = static_cast<GoldenPatternWithStat*>(omtfCand.getGoldenPatern());
    if( omtfCandGp->key().thePt > 100 && exptCandGp->key().thePt <= 15 ) {
      //edm::LogVerbatim("l1tOmtfEventPrint") <<iEvent.id()<<std::endl;
      cout<<" ptSim "<<ptSim<<" chargeSim "<<chargeSim<<std::endl;
      edm::LogVerbatim("l1tOmtfEventPrint") <<"iProcessor "<<iProcessor<<" exptCandGp "<<exptCandGp->key()<<std::endl;
      edm::LogVerbatim("l1tOmtfEventPrint") <<"iProcessor "<<iProcessor<<" omtfCandGp "<<omtfCandGp->key()<<std::endl;
      edm::LogVerbatim("l1tOmtfEventPrint") <<"omtfResult "<<std::endl<<omtfResult<<std::endl;
      int refHitNum = omtfCand.getRefHitNumber();
      edm::LogVerbatim("l1tOmtfEventPrint") <<"other gps results"<<endl;
      for(auto& gp : goldenPatterns) {
        if(omtfResult.getFiredLayerCnt() == gp->getResults()[procIndx][iRefHit].getFiredLayerCnt() )
        {
          edm::LogVerbatim("l1tOmtfEventPrint") <<gp->key()<<std::endl<<gp->getResults()[procIndx][iRefHit]<<std::endl;
        }
      }
      std::cout<<std::endl;
    }
  }*/
}

void EmulationObserverBase::observeEventBegin(const edm::Event& iEvent) {
  omtfCand.reset(new AlgoMuon());
  candProcIndx = 0xffff;
  //exptResult =  GoldenPatternResult();

  simMuon = findSimMuon(iEvent);
  //edm::LogVerbatim("l1tOmtfEventPrint") <<__FUNCTION__<<":"<<__LINE__<<" evevt "<<iEvent.id().event()<<" simMuon pt "<<simMuon->momentum().pt()<<" GeV "<<std::endl;
}

const SimTrack* EmulationObserverBase::findSimMuon(const edm::Event& event, const SimTrack* previous) {
  const SimTrack* result = nullptr;
  if (edmCfg.exists("simTracksTag") == false)
    return result;

  edm::Handle<edm::SimTrackContainer> simTks;
  event.getByLabel(edmCfg.getParameter<edm::InputTag>("simTracksTag"), simTks);

  //LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<" simTks->size() "<<simTks->size()<<std::endl;
  for (std::vector<SimTrack>::const_iterator it = simTks->begin(); it < simTks->end(); it++) {
    const SimTrack& aTrack = *it;
    if (!(aTrack.type() == 13 || aTrack.type() == -13))
      continue;
    if (previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(), previous->momentum()) < 0.07)
      continue;
    if (!result || aTrack.momentum().pt() > result->momentum().pt())
      result = &aTrack;
  }
  return result;
}

const std::vector<const reco::GenParticle*> EmulationObserverBase::findGenMuon(const edm::Event &event) {
  std::vector<const reco::GenParticle*> muons;

  if(edmCfg.exists("genParticleTag") == false)
    return muons;

  edm::Handle<reco::GenParticleCollection> genParticles;

  event.getByLabel(edmCfg.getParameter<edm::InputTag>("genParticleTag"), genParticles);

  //todo
  float etaCutFrom = 0.8;
  float etaCutTo = 1.24;

  for(size_t i = 0; i < genParticles->size() ; ++ i) {
    const reco::GenParticle& genPart = (*genParticles)[i];

    if( abs(genPart.pdgId() )== 13 ) {
      if( abs (genPart.momentum().eta() ) >= etaCutFrom && abs (genPart.momentum().eta() ) <= etaCutTo) {
        genPart.momentum().eta();

        muons.push_back( &genPart);
        //int id = p.pdgId();
        LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
                  <<" px "<<genPart.px()<<" py "<<genPart.py()<<" pz "<<genPart.pz()<<" pt "<<genPart.pt()
                  <<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<" eta "<<genPart.momentum().eta()<<endl;

        LogTrace("l1tOmtfEventPrint")<<__FUNCTION__<<":"<<__LINE__<<" genParticle: pdgId "<<genPart.pdgId()
                  <<" px "<<genPart.px()/genPart.p()<<" py "<<genPart.py()/genPart.p()<<" pz "<<genPart.pz()/genPart.p()<<" pt "<<genPart.pt()<<endl;
        //<<" vx "<<genPart.vx()<<" vy "<<genPart.vy()<<" vz "<<genPart.vz()<<endl;
      }
    }
  }

  return muons;
}
