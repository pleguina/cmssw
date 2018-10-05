/*
 * PdfGeneratorProcessor.cpp
 *
 *  Created on: Jul 28, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PdfGeneratorProcessor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdfGen.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
template<class GoldenPatternType>
void PdfGeneratorProcessor<GoldenPatternType>::fillCounts(unsigned int iProcessor,
    const OMTFinput & aInput,
    const SimTrack* aSimMuon) {
  std::ostringstream myStr;
  int theCharge = (abs(aSimMuon->type()) == 13) ? aSimMuon->type()/-13 : 0;
  myStr<<"aSimMuon->momentum().pt() "<<aSimMuon->momentum().pt()<<" theCharge "<<theCharge<<std::endl;
/*  unsigned int  iPt =  RPCConst::iptFromPt(aSimMuon->momentum().pt());
  myStr<<"line "<<__LINE__<<" iPt "<<iPt<<" theCharge "<<theCharge<<std::endl;
  ///Stupid conersion. Have to go through PAC pt scale, as we later
  ///shift resulting pt code by +1
  iPt+=1;
  if(iPt>31) iPt=200*2+1;
  else iPt = RPCConst::ptFromIpt(iPt)*2.0+1;//MicroGMT has 0.5 GeV step size, with lower bin edge  (uGMT_pt_code - 1)*step_size
  myStr<<"line "<<__LINE__<<" iPt "<<iPt<<std::endl;*/
  //////

  //////////////////////////////////////
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if(refHitsBits.none())
    return;

  myStr<<"iProcessor: "<<iProcessor<<std::endl;
  myStr<<"Input: ------------"<<std::endl;
  myStr<<aInput<<std::endl;
  edm::LogInfo("OMTF processor")<<myStr.str();

  //std::cout<<__FUNCTION__<<":"<<__LINE__<<" muon iPt "<<aSimMuon->momentum().pt()<<" theCharge "<<theCharge<<std::endl;
  for(unsigned int iLayer=0;iLayer< this->myOmtfConfig->nLayers();++iLayer) {

    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);
    if(!layerHits.size())
      continue;

    ///Number of reference hits to be checked.
    for(unsigned int iRefHit=0;iRefHit < this->myOmtfConfig->nRefHits();++iRefHit){

      if(!refHitsBits[iRefHit]) continue;
      const RefHitDef & aRefHitDef = this->myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int refLayerLogicNumber = this->myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      int phiRef = aInput.getLayerData(refLayerLogicNumber)[aRefHitDef.iInput];

      int refLayerPhiB = 0;
      if(refLayerLogicNumber < 6) {//is DT layer TODO - check
        refLayerPhiB = aInput.getLayerData(refLayerLogicNumber+1)[aRefHitDef.iInput]; //corresponding bending layer has number +1 versus the phi DT layer - this is configured in the XML
        if(this->myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) == 0) {
          throw cms::Exception("not good: the layer is not bending layer");
        }
      }


      unsigned int iRegion = aRefHitDef.iRegion;
      if(this->myOmtfConfig->getBendingLayers().count(iLayer))
        phiRef = 0;
      const OMTFinput::vector1D restrictedLayerHits = this->restrictInput(iProcessor, iRegion, iLayer,layerHits);
      for(auto& itGP: this->theGPs) {
/*        if(itGP->key().theCharge != theCharge)  //TODO it was in the orginal code, i commented it out, check why
          continue;*/

        if(aSimMuon->momentum().pt() >= this->myOmtfConfig->getPatternPtRange(itGP->key().number()).ptFrom &&
           aSimMuon->momentum().pt() < this->myOmtfConfig->getPatternPtRange(itGP->key().number()).ptTo    &&
           itGP->key().theCharge == theCharge ) {
          itGP->addCount(aRefHitDef.iRefLayer, iLayer, phiRef, restrictedLayerHits, refLayerPhiB);
          break;
        }
        //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iLayer "<<iLayer<<" refLayerPhiB "<<refLayerPhiB<<std::endl;
      }
    }
  }
}

///////////////////////////////////////////////
///////////////////////////////////////////////
template<class GoldenPatternType>
void  PdfGeneratorProcessor<GoldenPatternType>::averagePatterns(){
/*  OMTFConfiguration::vector2D mergedPartters = this->myOmtfConfig->getMergedPartters();
  for(unsigned int iGroup = 0; iGroup < mergedPartters.size(); iGroup++) {
    if(mergedPartters[iGroup].size() > 1) {

      GoldenPatternType* gp = this->theGPs.at(mergedPartters[iGroup][0]).get();
      if(gp == 0)
        throw cms::Exception("PdfGeneratorProcessor::averagePatterns works only for GoldenPatter");

      GoldenPattern::meanDistPhiArrayType meanDistPhi = gp->getMeanDistPhi();

      for(unsigned int i = 1; i < mergedPartters[iGroup].size(); i++) {
        GoldenPatternType* gp = this->theGPs.at( mergedPartters[iGroup][i]).get();
        for(unsigned int iLayer=0; iLayer < this->myOmtfConfig->nLayers(); ++iLayer) {
          for(unsigned int iRefLayer=0; iRefLayer < this->myOmtfConfig->nRefLayers(); ++iRefLayer) {
            meanDistPhi[iLayer][iRefLayer][0] += gp->getMeanDistPhi()[iLayer][iRefLayer][0];
          }
        }
      }

      for(unsigned int iLayer=0; iLayer< this->myOmtfConfig->nLayers(); ++iLayer) {
        for(unsigned int iRefLayer=0; iRefLayer < this->myOmtfConfig->nRefLayers(); ++iRefLayer) {
          meanDistPhi[iLayer][iRefLayer][0] /= mergedPartters[iGroup].size();
        }
      }

      for(unsigned int i = 0; i < mergedPartters[iGroup].size(); i++) {
        GoldenPatternType* gp = this->theGPs.at( mergedPartters[iGroup][i]).get();
        //shiftGP(gp, meanDistPhi, gp->getMeanDistPhi()); TODO FIXME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        gp->setMeanDistPhi(meanDistPhi);
      }
    }
  }*/
}
///////////////////////////////////////////////
///////////////////////////////////////////////
/*template<class GoldenPatternType>
void PdfGeneratorProcessor<GoldenPatternType>::shiftGP(GoldenPattern *aGP,
    const GoldenPattern::meanDistPhiArrayType& meanDistPhiNew,
    const GoldenPattern::meanDistPhiArrayType& meanDistPhiOld){

  ///Shift pdfs by differecne between original menaDistPhi, and
  ///the averaged value
  unsigned int nPdfBins =  exp2(myOmtfConfig->nPdfAddrBits());
  GoldenPattern::pdfArrayType pdfAllRef = aGP->getPdf();

  int indexShift = 0;
  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
    for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
      indexShift = meanDistPhiOld[iLayer][iRefLayer][0] - meanDistPhiNew[iLayer][iRefLayer][0];
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin)
        pdfAllRef[iLayer][iRefLayer][iPdfBin] = 0;
      for(unsigned int iPdfBin=0;iPdfBin<nPdfBins;++iPdfBin){
        if((int)(iPdfBin)+indexShift>=0 && iPdfBin+indexShift<nPdfBins)
          pdfAllRef[iLayer][iRefLayer][iPdfBin+indexShift] = aGP->pdfValue(iLayer, iRefLayer, iPdfBin);
      }
    }
  }
  aGP->setPdf(pdfAllRef);
}*/

template class PdfGeneratorProcessor<GoldenPatternPdf4D>;
template class PdfGeneratorProcessor<GoldenPatternPdfGen>;
