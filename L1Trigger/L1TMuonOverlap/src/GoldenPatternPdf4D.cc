#include <iostream>
#include <iomanip>
#include <cmath>

#include <TF1.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"


////////////////////////////////////////////////////
////////////////////////////////////////////////////
GoldenPatternPdf4D::GoldenPatternPdf4D(const Key & aKey, const OMTFConfiguration * omtfConfig) : IGoldenPattern(aKey, omtfConfig),
  theKey(aKey), myOmtfConfig(omtfConfig) {
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
IGoldenPattern::layerResult GoldenPatternPdf4D::process1Layer1RefLayer(unsigned int iRefLayer,
    unsigned int iLayer,
    const int phiRefHit,
    const OMTFinput::vector1D & layerHits,
    int refLayerPhiB){

  IGoldenPattern::layerResult aResult;

  int phiMean = meanDistPhi[iLayer][iRefLayer];
  int phiDist = exp2(myOmtfConfig->nPdfAddrBits());
  ///Select hit closest to the mean of probability 
  ///distribution in given layer
  for(auto itHit: layerHits){
    if(itHit>=(int)myOmtfConfig->nPhiBins()) continue;    
    if(abs(itHit-phiMean-phiRefHit) < abs(phiDist))
      phiDist = itHit-phiMean-phiRefHit;
  }

  ///Check if phiDist is within pdf range -63 +63 
  if(abs(phiDist) > (exp2(myOmtfConfig->nPdfAddrBits()-1) -1))
    return aResult;

  ///Shift phidist, so 0 is at the middle of the range
  phiDist+=exp2(myOmtfConfig->nPdfAddrBits()-1);

  int pdfVal = pdfAllRef[iLayer][iRefLayer][refLayerPhiB][phiDist];

  return IGoldenPattern::layerResult(pdfVal,pdfVal>0);
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
int GoldenPatternPdf4D::propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer){

  unsigned int iLayer = 2; //MB2 
  //if(etaRef>101) iLayer = 7;//RE2
  return phiRef + meanDistPhi[iLayer][iRefLayer];

}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternPdf4D::addCount(unsigned int iRefLayer,
    unsigned int iLayer,
    const int phiRefHit,
    const OMTFinput::vector1D & layerHits,
    int refLayerPhiB) {

  int nHitsInLayer = 0;
  int phiDist = 1<<myOmtfConfig->nPdfAddrBits();//=128
  for(auto itHit: layerHits) {
    if(itHit >= (int)myOmtfConfig->nPhiBins())
      continue; //why it can ever happened???
    if(abs(itHit-phiRefHit) < phiDist)
      phiDist = itHit-phiRefHit;
    ++nHitsInLayer;
  }
  ///For making the patterns take events with a single hit in each layer
  if(nHitsInLayer>1 || nHitsInLayer==0) return;

  if(pdfAllRef[iLayer][iRefLayer].size()  > 1) //the ref layer is bending layer, so the it is possible that refLayerPhiB != 0
    refLayerPhiB = refLayerPhiB + (1<<(myOmtfConfig->nPdfAddrBits()-1));
  if(refLayerPhiB < 0 || (unsigned int)refLayerPhiB >= pdfAllRef[iLayer][iRefLayer].size() ) {
    edm::LogError("GoldenPatternPdf4D::addCount ")<<"iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" refLayerPhiB "<<refLayerPhiB
        <<" out of bounds. pdfAllRef[iLayer][iRefLayer].size() "<<pdfAllRef[iLayer][iRefLayer].size();
    return;
  }

  ///Shift phiDist so it is in +-Pi range
  if(phiDist>=(int)myOmtfConfig->nPhiBins()/2) phiDist-=(int)myOmtfConfig->nPhiBins();
  if(phiDist<=-(int)myOmtfConfig->nPhiBins()/2) phiDist+=(int)myOmtfConfig->nPhiBins();

  ///Shift phidist, so 0 is at the middle of the range
  int phiDistShift=phiDist+exp2(myOmtfConfig->nPdfAddrBits()-1);

  ///Check if phiDist is within pdf range
  ///in -64 +63 U2 code
  ///Find more elegant way to check this.
  if(phiDistShift < 0 || (unsigned int)phiDistShift >= pdfAllRef[iLayer][iRefLayer][refLayerPhiB].size() ) {
    edm::LogError("GoldenPatternPdf4D::addCount ")<<"phiDistShift "<<phiDistShift<<" out of bounds";
    return;
  }


  if((int)iLayer == myOmtfConfig->getRefToLogicNumber()[iRefLayer])
    ++meanDistPhiCounts[iLayer][iRefLayer];

  pdfAllRef[iLayer][iRefLayer].at(refLayerPhiB).at(phiDistShift)++;

/*  if(linearFitters.at(iLayer).at(iRefLayer) != 0 ) {
    double refLayerPhiB_d = refLayerPhiB;
    linearFitters.at(iLayer).at(iRefLayer)->AddPoint(&refLayerPhiB_d, phiDistShift);
    //std::cout<<__FUNCTION__<<" "<<refLayerPhiB_d<<" "<<phiDistShift<<std::endl;
  }*/
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const GoldenPatternPdf4D & aPattern){

  out <<"GoldenPatternPdf4D "<< aPattern.theKey <<std::endl;
  out <<"Number of reference layers: "<<aPattern.meanDistPhi[0].size()
          <<", number of measurement layers: "<<aPattern.pdfAllRef.size()
          <<std::endl;

  if(!aPattern.meanDistPhi.size()) return out;
  if(!aPattern.pdfAllRef.size()) return out;

  out<<"Mean dist phi per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer<<" (";
    for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
      out<<std::setw(3)<<aPattern.meanDistPhi[iLayer][iRefLayer]<<"\t";
    }
    out<<")"<<std::endl;
  }

  if(aPattern.meanDistPhiCounts.size()){
    out<<"Counts number per layer:"<<std::endl;
    for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
      out<<"Ref layer: "<<iRefLayer<<" (";
      for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
        out<<aPattern.meanDistPhiCounts[iLayer][iRefLayer]<<"\t";
      }
      out<<")"<<std::endl;
    }
  }

  unsigned int nPdfAddrBits = 7;
  out<<"PDF per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.pdfAllRef[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer;
    for (unsigned int iLayer=0;iLayer<aPattern.pdfAllRef.size();++iLayer){   
      out<<", measurement layer: "<<iLayer<<std::endl;
      for (unsigned int iRefPhiB=0; iRefPhiB < aPattern.pdfAllRef[iLayer][iRefLayer].size(); ++iRefPhiB) {
        for(unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){
          out<<std::setw(2)<<aPattern.pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf]<<" ";
        }
      }
      out<<std::endl;
    }
  }

  return out;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternPdf4D::reset(){

  IGoldenPattern::vector1D meanDistPhi1D(myOmtfConfig->nRefLayers());
  IGoldenPattern::vector2D meanDistPhi2D(myOmtfConfig->nLayers());
  meanDistPhi2D.assign(myOmtfConfig->nLayers(), meanDistPhi1D);
  meanDistPhi = meanDistPhi2D;
  meanDistPhiCounts = meanDistPhi2D;

  unsigned int nBins = exp2(myOmtfConfig->nPdfAddrBits());
  IGoldenPattern::vector1D pdf1D(nBins);
  //IGoldenPattern::vector2D pdf2D(myOmtfConfig->nRefLayers());
  IGoldenPattern::vector3D pdf3D(myOmtfConfig->nRefLayers());
  //GoldenPatternPdf4D::vector4D pdf4D(myOmtfConfig->nLayers());

  std::cout<<"GoldenPatternPdf4D::reset() myOmtfConfig->nPdfAddrBits() "<<myOmtfConfig->nPdfAddrBits()<<" nBins "<<nBins<<std::endl;
  std::cout<<"GoldenPatternPdf4D::reset() myOmtfConfig->nLayers() "<<myOmtfConfig->nLayers()<<" myOmtfConfig->nRefLayers()  "<<myOmtfConfig->nRefLayers()<<std::endl;

  linearFitters.assign(myOmtfConfig->nLayers(), std::vector<TLinearFitter*>(myOmtfConfig->nRefLayers(), 0) );
  linearFits.assign(myOmtfConfig->nLayers(), std::vector<TF1*>(myOmtfConfig->nRefLayers(), 0) );;

  pdfAllRef.assign(myOmtfConfig->nLayers(), pdf3D);
  for (unsigned int iLayer=0; iLayer<pdfAllRef.size(); ++iLayer) {
    for (unsigned int iRefLayer=0; iRefLayer<pdfAllRef[0].size(); ++iRefLayer) {
      int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[iRefLayer];
      if(myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) != 0) {//ref leyer is DT layer with banding
        pdfAllRef[iLayer][iRefLayer].assign(nBins, pdf1D);

        linearFitters[iLayer][iRefLayer] = new TLinearFitter();
        linearFitters[iLayer][iRefLayer]->StoreData(kTRUE);
        std::ostringstream ostr;
        ostr<<"fit_lay_"<<iLayer<<"_refLay_"<<iRefLayer;
        linearFits[iLayer][iRefLayer] = new TF1(ostr.str().c_str(), "pol1", 0, nBins);
        linearFitters[iLayer][iRefLayer]->SetFormula(linearFits[iLayer][iRefLayer]->GetFormula());
      }
      else {
        pdfAllRef[iLayer][iRefLayer].assign(1, pdf1D);
      }
      std::cout<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" refLayerLogicNumber "<<refLayerLogicNumber<<" pdfAllRef[iLayer][iRefLayer].size() "
          <<pdfAllRef[iLayer][iRefLayer].size()<<" pdfAllRef[iLayer][iRefLayer][0].size() "<<pdfAllRef[iLayer][iRefLayer][0].size()<<std::endl;
    }
  }
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPatternPdf4D::normalise(unsigned int nPdfAddrBits) {
    std::cout<<__FUNCTION__<<" "<<key()<<std::endl;
/*    for(unsigned int iLayer = 0; iLayer < linearFitters.size(); iLayer++) {
      for(unsigned int iRefLayer=0; iRefLayer < linearFitters.at(iLayer).size(); iRefLayer++) {
        if(linearFitters.at(iLayer).at(iRefLayer) != 0 ) {
          int res = linearFitters.at(iLayer).at(iRefLayer)->Eval();
          std::cout<<"iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" linearFitter result "<<res;
          if(res == 0) {
            std::cout<<linearFits[iLayer][iRefLayer]->GetParName(0)<<" "<<linearFits[iLayer][iRefLayer]->GetParameter(0)<<" ";
            std::cout<<linearFits[iLayer][iRefLayer]->GetParName(1)<<" "<<linearFits[iLayer][iRefLayer]->GetParameter(1);
          }
          std::cout<<std::endl;
        }
      }
    }*/

    for(unsigned int iLayer = 0; iLayer < linearFitters.size(); iLayer++) {
      for(unsigned int iRefLayer=0; iRefLayer < linearFitters.at(iLayer).size(); iRefLayer++) {
        if(linearFitters.at(iLayer).at(iRefLayer) != 0 ) {
          for(unsigned int iRefLayerPhiB = 0; iRefLayerPhiB < pdfAllRef[iLayer][iRefLayer].size(); iRefLayerPhiB++) {
            double sum_p = 0;
            double sum_p_x = 0;
            for(unsigned int iLayerPhi=0; iLayerPhi < pdfAllRef[iLayer][iRefLayer][iRefLayerPhiB].size(); iLayerPhi++) {
              sum_p += pdfAllRef[iLayer][iRefLayer][iRefLayerPhiB][iLayerPhi];
              sum_p_x += pdfAllRef[iLayer][iRefLayer][iRefLayerPhiB][iLayerPhi] * iLayerPhi;
            }

            if(sum_p > 50) {
              double meanPhi = sum_p_x/sum_p;
              double x = iRefLayerPhiB;
              linearFitters.at(iLayer).at(iRefLayer)->AddPoint(&x, meanPhi);
            }
          }

          int res = linearFitters.at(iLayer).at(iRefLayer)->Eval();
          std::cout<<"iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" linearFitter result "<<res;
          if(res == 0) {
            std::cout<<linearFits[iLayer][iRefLayer]->GetParName(0)<<" "<<linearFits[iLayer][iRefLayer]->GetParameter(0)<<" ";
            std::cout<<linearFits[iLayer][iRefLayer]->GetParName(1)<<" "<<linearFits[iLayer][iRefLayer]->GetParameter(1);
          }
          std::cout<<std::endl;
        }
      }
    }
///not needed here
/*  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer) {
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer) {
      for (unsigned int iRefPhiB=0; iRefPhiB < pdfAllRef[iLayer][iRefLayer].size(); ++iRefPhiB) {
        for (unsigned int iPdf=0;iPdf<pdfAllRef[iLayer][iRefLayer][iRefPhiB].size();++iPdf) {
          float pVal = log((float)pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf] / meanDistPhiCounts[myOmtfConfig->getRefToLogicNumber()[iRefLayer]][iRefLayer]);
          if(pVal < log(myOmtfConfig->minPdfVal()))
            continue;
          meanDistPhi[iLayer][iRefLayer]+=(iPdf - exp2(myOmtfConfig->nPdfAddrBits()-1))*pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf];
          if((int)iLayer!=myOmtfConfig->getRefToLogicNumber()[iRefLayer]) //this iLayer is not the refLayer
            meanDistPhiCounts[iLayer][iRefLayer]+=pdfAllRef[iLayer][iRefLayer][iRefPhiB][iPdf];
        }
      }
    }
  }

  ///Mean dist phi  
  for (unsigned int iRefLayer=0;iRefLayer<meanDistPhi[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<meanDistPhi.size();++iLayer){   
      if(meanDistPhiCounts.size() && meanDistPhiCounts[iLayer][iRefLayer]){
        if(meanDistPhiCounts[iLayer][iRefLayer]<1000)
          meanDistPhi[iLayer][iRefLayer] = 0;
        else
          meanDistPhi[iLayer][iRefLayer] = rint((float)meanDistPhi[iLayer][iRefLayer]/meanDistPhiCounts[iLayer][iRefLayer]);
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

  vector3D  pdfAllRefTmp = pdfAllRef;
  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer){   
      for (unsigned int iPdf=0;iPdf<pdfAllRef[iLayer][iRefLayer].size();++iPdf){
        pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
        ///Shift pdf index by meanDistPhi
        int index = iPdf - exp2(myOmtfConfig->nPdfAddrBits()-1)  - meanDistPhi[iLayer][iRefLayer] + exp2(nPdfAddrBits-1);
        if(index<0 || index>exp2(nPdfAddrBits)-1) continue;
        pdfAllRef[iLayer][iRefLayer][index] = pdfAllRefTmp[iLayer][iRefLayer][iPdf];
      }
    }
  }*/
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
bool GoldenPatternPdf4D::hasCounts(){

  for (unsigned int iRefLayer=0;iRefLayer<meanDistPhi[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<meanDistPhi.size();++iLayer){   
      if(meanDistPhiCounts.size() && meanDistPhiCounts[iLayer][iRefLayer]) return true;
    }
  }
  return false;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
