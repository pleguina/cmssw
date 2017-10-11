/*
 * ProcessorBase.hxx
 *
 *  Created on: Aug 1, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PROCESSORBASE_HXX_
#define OMTF_PROCESSORBASE_HXX_

#include <memory>

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

template<class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::resetConfiguration(){
  //myResults.clear();
  //for(auto it: theGPs) delete it;
  theGPs.clear();
}

///////////////////////////////////////////////
///////////////////////////////////////////////
template<class GoldenPatternType>
bool ProcessorBase<GoldenPatternType>::configure(const OMTFConfiguration* omtfConfig,
    const L1TMuonOverlapParams* omtfPatterns){
  resetConfiguration();

  myOmtfConfig = omtfConfig;

  //myResults.assign(myOmtfConfig->nTestRefHits(),ProcessorBase::resultsMap());

  const l1t::LUT* chargeLUT =  omtfPatterns->chargeLUT();
  const l1t::LUT* etaLUT =  omtfPatterns->etaLUT();
  const l1t::LUT* ptLUT =  omtfPatterns->ptLUT();
  const l1t::LUT* pdfLUT =  omtfPatterns->pdfLUT();
  const l1t::LUT* meanDistPhiLUT =  omtfPatterns->meanDistPhiLUT();

  unsigned int nGPs = myOmtfConfig->nGoldenPatterns();
  edm::LogInfo("MTFProcessor::configure")<<"myOmtfConfig->nGoldenPatterns() "<<nGPs<<std::endl;
  unsigned int address = 0;
  unsigned int iEta, iPt;
  int iCharge;
  int meanDistPhiSize = myOmtfConfig->nLayers() * myOmtfConfig->nRefLayers() * myOmtfConfig->nGoldenPatterns();
  for(unsigned int iGP=0;iGP<nGPs;++iGP){
    address = iGP;
    iEta = etaLUT->data(address);
    iCharge = chargeLUT->data(address)==0? -1:1;
    iPt = ptLUT->data(address);

    Key aKey(iEta,iPt,iCharge,iGP);
    edm::LogInfo("ProcessorBase::configure")<<"adding pattern "<<aKey<<" "<<myOmtfConfig->getPatternPtRange(iGP).ptFrom<<" - "<<myOmtfConfig->getPatternPtRange(iGP).ptTo<<" GeV"<<std::endl;
    GoldenPatternType* aGP = new GoldenPatternType(aKey, myOmtfConfig);

    ///Mean dist phi data
    for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        address = iRefLayer + iLayer*myOmtfConfig->nRefLayers() + iGP*(myOmtfConfig->nRefLayers()*myOmtfConfig->nLayers());
        int value = meanDistPhiLUT->data(address) - (1<<(meanDistPhiLUT->nrBitsData() -1));
        aGP->setMeanDistPhiValue(value, iLayer, iRefLayer, 0);
        if(meanDistPhiLUT->nrBitsAddress() == 15) {//for the new version of the meanDistPhi which have two values for each gp,iLayer,iRefLayer, FIXME: do it a better way
          value = meanDistPhiLUT->data(address + meanDistPhiSize) - (1<<(meanDistPhiLUT->nrBitsData() -1));
          //the second meanDistPhi is in the LUT at the position (address+meanDistPhiSize)
          aGP->setMeanDistPhiValue(value, iLayer, iRefLayer, 1);
        }
      }
      ///Pdf data
      for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer){
        for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf){
          address = iPdf + iRefLayer*(1<<myOmtfConfig->nPdfAddrBits()) +
              iLayer*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits()) +
              iGP*myOmtfConfig->nLayers()*myOmtfConfig->nRefLayers()*(1<<myOmtfConfig->nPdfAddrBits());
          int value = pdfLUT->data(address);
          aGP->setPdfValue(value, iLayer, iRefLayer, iPdf);
        }
      }
    }
    addGP(aGP);
  }

  return true;
}

///////////////////////////////////////////////
///////////////////////////////////////////////
template<class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::addGP(GoldenPatternType* aGP) {
  theGPs.emplace_back(std::unique_ptr<GoldenPatternType>(aGP));
}

////////////////////////////////////////////
////////////////////////////////////////////
template<class GoldenPatternType>
OMTFinput::vector1D ProcessorBase<GoldenPatternType>::restrictInput(unsigned int iProcessor,
    unsigned int iRegion,
    unsigned int iLayer,
    const OMTFinput::vector1D & layerHits) {

  OMTFinput::vector1D myHits = layerHits;

  unsigned int iStart = myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].first;
  unsigned int iEnd = iStart + myOmtfConfig->getConnections()[iProcessor][iRegion][iLayer].second -1;

  for(unsigned int iInput=0;iInput<myHits.size();++iInput){
    if(iInput<iStart || iInput>iEnd) myHits[iInput] = myOmtfConfig->nPhiBins();
  }
  return myHits;
}

////////////////////////////////////////////
////////////////////////////////////////////
template<class GoldenPatternType>
void ProcessorBase<GoldenPatternType>::fillCounts(unsigned int iProcessor,
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
  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer) {

    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);
    if(!layerHits.size())
      continue;

    ///Number of reference hits to be checked.
    for(unsigned int iRefHit=0;iRefHit<myOmtfConfig->nRefHits();++iRefHit){

      if(!refHitsBits[iRefHit]) continue;
      const RefHitDef & aRefHitDef = myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      int phiRef = aInput.getLayerData(refLayerLogicNumber)[aRefHitDef.iInput];

      int refLayerPhiB = 0;
      if(refLayerLogicNumber < 6) {//is DT layer TODO - check
        refLayerPhiB = aInput.getLayerData(refLayerLogicNumber+1)[aRefHitDef.iInput]; //corresponding bending layer has number +1 versus the phi DT layer - this is configured in the XML
        if(myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) == 0) {
          throw cms::Exception("not good: the layer is not bending layer");
        }
      }


      unsigned int iRegion = aRefHitDef.iRegion;
      if(myOmtfConfig->getBendingLayers().count(iLayer))
        phiRef = 0;
      const OMTFinput::vector1D restrictedLayerHits = restrictInput(iProcessor, iRegion, iLayer,layerHits);
      for(auto& itGP: theGPs) {
/*        if(itGP->key().theCharge != theCharge)  //TODO it was in the orginal code, i commented it out, check why
          continue;*/

        if(aSimMuon->momentum().pt() >= myOmtfConfig->getPatternPtRange(itGP->key().number()).ptFrom &&
           aSimMuon->momentum().pt() < myOmtfConfig->getPatternPtRange(itGP->key().number()).ptTo    &&
           itGP->key().theCharge == theCharge ) {
          itGP->addCount(aRefHitDef.iRefLayer, iLayer, phiRef, restrictedLayerHits, refLayerPhiB);
          break;
        }
        //std::cout<<__FUNCTION__<<":"<<__LINE__<<" iLayer "<<iLayer<<" refLayerPhiB "<<refLayerPhiB<<std::endl;
      }
    }
  }
}

#endif /* OMTF_PROCESSORBASE_HXX_ */
