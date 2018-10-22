/*
 * PatternGeneratorTT.cc
 *
 *  Created on: Oct 17, 2018
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PatternGeneratorTT.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessorTTMerger.h"

#include "TFile.h"

PatternGeneratorTT::PatternGeneratorTT(const edm::ParameterSet& edmCfg, const OMTFConfiguration* omtfConfig, std::vector<std::shared_ptr<GoldenPatternWithStat> >& gps):
  PatternOptimizerBase(edmCfg, omtfConfig, gps)
{
  // TODO Auto-generated constructor stub

  for(auto& gp : goldenPatterns) {
      if(gp->key().thePt == 0)
        continue;

      gp->iniStatisitics(gp->getPdf()[0][0].size() * 4, 1); //TODO
  }


  //TODO remove, it is needed only for the patterns generation
  for(auto& gp : goldenPatterns) {
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
        if(refLayerLogicNum == iLayer) {
          for(unsigned int iBin = 0; iBin < gp->getPdf()[iLayer][iRefLayer].size(); iBin++) {
            gp->pdfAllRef[iLayer][iRefLayer][iBin] = 1;
          }
        }
      }
    }
  }

  ptDeltaPhiHists.resize(2);
  for(unsigned int iCharge = 0; iCharge <= 1; iCharge++) {
    for(unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) { //for the moment filing only ref layer, remove whe
      if(iLayer == 0 || iLayer == 2 || iLayer == 4 || iLayer == 6 || iLayer == 7 || iLayer == 10 || iLayer == 11  || iLayer == 16) {
        ostringstream name;
        name<<"ptDeltaPhiHist_ch_"<<iCharge<<"_Layer_"<<iLayer;
        int phiFrom = -10;
        int phiTo   = 1000;
        int phiBins = phiTo - phiFrom;

        if(iCharge == 1) {
          phiFrom = -1000;
          phiTo = 10;
        }

        TH2I* ptDeltaPhiHist = new TH2I(name.str().c_str(), name.str().c_str(), 400, 0, 200, phiBins, phiFrom -0.5, phiTo -0.5);
        //cout<<"BinLowEdge "<<ptDeltaPhiHist->GetYaxis()->GetBinLowEdge(100)<<" BinUpEdge "<<ptDeltaPhiHist->GetYaxis()->GetBinUpEdge(100);
        ptDeltaPhiHists[iCharge].push_back(ptDeltaPhiHist);
      }
      else
        ptDeltaPhiHists[iCharge].push_back(nullptr);
    }
  }
}

PatternGeneratorTT::~PatternGeneratorTT() {
  // TODO Auto-generated destructor stub
}

void PatternGeneratorTT::updateStat() {
  cout<<__FUNCTION__<<":"<<__LINE__<<" omtfCand "<<*omtfCand<<std::endl;;
  TTAlgoMuon* ttAlgoMuon = dynamic_cast<TTAlgoMuon*>(omtfCand.get());
  if(!ttAlgoMuon) {
    cout<<__FUNCTION__<<":"<<__LINE__<<" ttAlgoMuon is null"<<std::endl;
    throw runtime_error("ttAlgoMuon is null");
  }
  //TODO go through all refLayers in the omtfCandGp

  GoldenPatternWithStat* omtfCandGp = static_cast<GoldenPatternWithStat*>(omtfCand->getGoldenPatern());

  unsigned int iCharge = omtfCand->getCharge();
  if(iCharge != 1)
    iCharge = 0;

  int pdfMiddle = 1<<(omtfConfig->nPdfAddrBits()-1);

  for(auto& gpResult : ttAlgoMuon->getGpResults()) {
    unsigned int iRefLayer =  gpResult.getRefLayer();
    unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];

    //if(gpResult.isValid() )
    {
      for(unsigned int iLayer = 0;  iLayer < gpResult.getHitPdfBins().size(); iLayer++) {
        //updating statistic for the gp which should have fired

        int phiDist = gpResult.getHitPdfBins()[iLayer];
        phiDist += omtfCandGp->meanDistPhiValue(iLayer, iRefLayer) - pdfMiddle;

        if(iLayer == refLayerLogicNum) //TODO remove if needed to far all layers
          ptDeltaPhiHists[iCharge][refLayerLogicNum]->Fill(ttAlgoMuon->getTtTrack().getPt(), phiDist);

        if( abs(phiDist) < omtfCandGp->getStatistics()[iLayer][iRefLayer].size() / 2 ) {
          phiDist += omtfCandGp->getStatistics()[iLayer][iRefLayer].size() / 2; //shifting so that phiDist = 0 is in the middle of the statisitics[iLayer][iRefLayer][]

          //updating statistic for the gp which found the candidate
          //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for omtf gp "<<omtfCandGp<<std::endl;
          omtfCandGp->updateStat(iLayer, iRefLayer, phiDist, 0, 1);
          //cout<<__FUNCTION__<<":"<<__LINE__<<" updating statistic for omtfCandGp: iLayer "<<iLayer<<" RefLayer "<<omtfCand.getRefLayer()<<" iBinOmtf "<<gpResult.getHitPdfBins()[iLayer]<<" fired "<<gpResult.isLayerFired(iLayer)<<std::endl;
        }

      }
    }
  }
}

void PatternGeneratorTT::observeEventEnd(const edm::Event& iEvent) {
  if(simMuon == 0 || omtfCand->getGoldenPatern() == 0)//no sim muon or empty candidate
    return;

  PatternOptimizerBase::observeEventEnd(iEvent);

  updateStat();
}

void PatternGeneratorTT::endJob() {

  upadatePdfRefLay();

  PatternOptimizerBase::endJob();
}

void PatternGeneratorTT::upadatePdfRefLay() {
  cout<<__FUNCTION__<<": "<<__LINE__<<endl;
  for(auto& gp : goldenPatterns) {
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {

        for(unsigned int iRefLayer = 0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
          unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[iRefLayer];
          if(refLayerLogicNum == iLayer) {

            //calculate meanDistPhi
            double meanDistPhi = 0;
            double count = 0;
            for(unsigned int iBin = 0; iBin < gp->getStatistics()[iLayer][iRefLayer].size(); iBin++) {
              meanDistPhi +=  iBin * gp->getStatistics()[iLayer][iRefLayer][iBin][0];
              count       +=         gp->getStatistics()[iLayer][iRefLayer][iBin][0];
            }
            if(count != 0) {
              meanDistPhi /= count;
              meanDistPhi -= gp->getStatistics()[iLayer][iRefLayer].size() / 2;


              cout<<__FUNCTION__<<": "<<__LINE__<<gp->key()<<" iLayer "<<iLayer<<" count "<<count<<" meanDistPhi "<<meanDistPhi<<endl;
              gp->setMeanDistPhiValue(meanDistPhi, iLayer, iRefLayer);
            }
          }
        }
      }
    }
}

void PatternGeneratorTT::saveHists(TFile& outfile) {
  outfile.mkdir("ptDeltaPhiHists" )->cd();
  for(unsigned int iCharge = 0; iCharge <= 1; iCharge++) {
    for(unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) { //for the moment filing only ref layer, remove whe
      if(ptDeltaPhiHists[iCharge][iLayer]) {
        ptDeltaPhiHists[iCharge][iLayer]->Write();
      }
    }
  }
}

