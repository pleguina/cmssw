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
GoldenPatternPdf4D::GoldenPatternPdf4D(const Key& aKey, const OMTFConfiguration * omtfConfig) : IGoldenPattern(aKey, omtfConfig) {
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

  if(pdfAllRef[iLayer][iRefLayer].size()  > 1) //the ref layer is bending layer, so it is possible that refLayerPhiB != 0
    refLayerPhiB = refLayerPhiB + (1<<(myOmtfConfig->nPdfAddrBits()-1)); //to have only positive values
  if(refLayerPhiB < 0 || (unsigned int)refLayerPhiB >= pdfAllRef[iLayer][iRefLayer].size() ) {
    edm::LogError("GoldenPatternPdf4D::addCount ")<<"iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" refLayerPhiB "<<refLayerPhiB
        <<" out of bounds. pdfAllRef[iLayer][iRefLayer].size() "<<pdfAllRef[iLayer][iRefLayer].size();
    return;
  }

  ///Shift phiDist so it is in +-Pi range
  if(phiDist>=(int)myOmtfConfig->nPhiBins()/2) phiDist-=(int)myOmtfConfig->nPhiBins();
  if(phiDist<=-(int)myOmtfConfig->nPhiBins()/2) phiDist+=(int)myOmtfConfig->nPhiBins();

  ///Shift phidist, so 0 is at the middle of the range
  int phiDistShift=phiDist + (1<<(myOmtfConfig->nPdfAddrBits()-1));

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

  pdfAllRef.assign(myOmtfConfig->nLayers(), pdf3D);
  for (unsigned int iLayer=0; iLayer<pdfAllRef.size(); ++iLayer) {
    for (unsigned int iRefLayer=0; iRefLayer<pdfAllRef[0].size(); ++iRefLayer) {
      int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[iRefLayer];
      if(myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) != 0) {//ref leyer is DT layer with banding
        pdfAllRef[iLayer][iRefLayer].assign(nBins, pdf1D);
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
