/*
 * GoldenPatternPdfGen.cc
 *
 *  Created on: Oct 23, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdfGen.h"

GoldenPatternPdfGen::GoldenPatternPdfGen(const Key & aKey, const OMTFConfiguration* omtfConfig):
  GoldenPattern(aKey, omtfConfig),
  meanDistPhiCounts(boost::extents[omtfConfig->nLayers()][omtfConfig->nRefLayers()]) {
  // TODO Auto-generated constructor stub

}

GoldenPatternPdfGen::~GoldenPatternPdfGen() {
  // TODO Auto-generated destructor stub
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternPdfGen::addCount(unsigned int iRefLayer,
    unsigned int iLayer,
    const int phiRefHit,
    const OMTFinput::vector1D & layerHits,
    int refLayerPhiB){

  int nHitsInLayer = 0;
  int phiDist = exp2(myOmtfConfig->nPdfAddrBits());//=128
  for(auto itHit: layerHits){
    if(itHit>=(int)myOmtfConfig->nPhiBins()) continue; //why it can ever happened???
    if(abs(itHit-phiRefHit)<phiDist) phiDist = itHit-phiRefHit;
    ++nHitsInLayer;
  }
  ///For making the patterns take events with a single hit in each layer
  if(nHitsInLayer>1 || nHitsInLayer==0) return;

  ///Shift phiDist so it is in +-Pi range
  if(phiDist>=(int)myOmtfConfig->nPhiBins()/2) phiDist-=(int)myOmtfConfig->nPhiBins();
  if(phiDist<=-(int)myOmtfConfig->nPhiBins()/2) phiDist+=(int)myOmtfConfig->nPhiBins();

  ///Shift phidist, so 0 is at the middle of the range
  int phiDistShift=phiDist+exp2(myOmtfConfig->nPdfAddrBits()-1);

  ///Check if phiDist is within pdf range
  ///in -64 +63 U2 code
  ///Find more elegant way to check this.
  if(phiDistShift<0 ||
      phiDistShift>exp2(myOmtfConfig->nPdfAddrBits())-1){
    return;
  }

  if((int)iLayer==myOmtfConfig->getRefToLogicNumber()[iRefLayer] )
    ++meanDistPhiCounts[iLayer][iRefLayer];

  ++pdfAllRef[iLayer][iRefLayer][phiDistShift];
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternPdfGen::normalise(unsigned int nPdfAddrBits){
  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer){
      for (unsigned int iPdf=0;iPdf<pdfAllRef[iLayer][iRefLayer].size();++iPdf){
        float pVal = log((float)pdfAllRef[iLayer][iRefLayer][iPdf]/meanDistPhiCounts[myOmtfConfig->getRefToLogicNumber()[iRefLayer]][iRefLayer]);
        if(pVal<log(myOmtfConfig->minPdfVal()))
          continue;
        meanDistPhi[iLayer][iRefLayer][0] +=(iPdf - exp2(myOmtfConfig->nPdfAddrBits()-1))*pdfAllRef[iLayer][iRefLayer][iPdf];
        if((int)iLayer!=myOmtfConfig->getRefToLogicNumber()[iRefLayer]) //this iLayer is not the refLayer
          meanDistPhiCounts[iLayer][iRefLayer]+=pdfAllRef[iLayer][iRefLayer][iPdf];
      }
    }
  }

  ///Mean dist phi
  for (unsigned int iRefLayer=0;iRefLayer<meanDistPhi[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<meanDistPhi.size();++iLayer){
      if(meanDistPhiCounts.size() && meanDistPhiCounts[iLayer][iRefLayer]){
        if(meanDistPhiCounts[iLayer][iRefLayer]<1000)
          meanDistPhi[iLayer][iRefLayer][0] = 0;
        else
          meanDistPhi[iLayer][iRefLayer][0] = rint((float)meanDistPhi[iLayer][iRefLayer][0]/meanDistPhiCounts[iLayer][iRefLayer]);
      }
    }
  }
  const float minPlog =  log(myOmtfConfig->minPdfVal());
  const unsigned int nPdfValBits = myOmtfConfig->nPdfValBits();
  ///Probabilities. Normalise and change from float to integer values
  float pVal;
  int digitisedVal, truncatedValue;
  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer){
      for (unsigned int iPdf=0;iPdf<pdfAllRef[iLayer][iRefLayer].size();++iPdf){
        if(!meanDistPhiCounts[myOmtfConfig->getRefToLogicNumber()[iRefLayer]][iRefLayer] ||
            !pdfAllRef[iLayer][iRefLayer][iPdf])
          continue;
        pVal = log((float)pdfAllRef[iLayer][iRefLayer][iPdf]/meanDistPhiCounts[myOmtfConfig->getRefToLogicNumber()[iRefLayer]][iRefLayer]);
        //the normalization is by all muons that crossed the refLayer, and not by all muons that crossed the given iLayer.
        //This is ok, since in a given iLayer the geometrical acceptance might depend on the pt,
        ///If there are only a few counts in given measurement layer, set pdf value to 0
        if((pVal<minPlog || meanDistPhiCounts[iLayer][iRefLayer]<1000)){
          pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
          continue;
        }
        ///Digitisation
        ///Values remapped 0->std::pow(2,nPdfValBits)
        ///          minPlog->0
        digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
        ///Make sure digitised value is saved using nBitsVal bits
        truncatedValue  = 0 | (digitisedVal  & ((int)pow(2,nPdfValBits)-1)); //FIXME why it is needed at all? if digitisedVal > ((int)pow(2,nPdfValBits)-1) the result is bad
        pdfAllRef[iLayer][iRefLayer][iPdf] = truncatedValue;
      }
    }
  }

  pdfArrayType pdfAllRefTmp = pdfAllRef;
  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer){
      for (unsigned int iPdf=0;iPdf<pdfAllRef[iLayer][iRefLayer].size();++iPdf){
        pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
        ///Shift pdf index by meanDistPhi
        int index = iPdf - exp2(myOmtfConfig->nPdfAddrBits()-1)  - meanDistPhi[iLayer][iRefLayer][0] + exp2(nPdfAddrBits-1);
        if(index<0 || index>exp2(nPdfAddrBits)-1) continue;
        pdfAllRef[iLayer][iRefLayer][index] = pdfAllRefTmp[iLayer][iRefLayer][iPdf];
      }
    }
  }
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
bool GoldenPatternPdfGen::hasCounts(){
  for (unsigned int iRefLayer=0;iRefLayer<meanDistPhiCounts[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<meanDistPhiCounts.size();++iLayer){
      if(meanDistPhiCounts.size() && meanDistPhiCounts[iLayer][iRefLayer]) return true;
    }
  }
  return false;
}
