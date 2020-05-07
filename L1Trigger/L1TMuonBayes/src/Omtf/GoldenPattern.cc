#include <L1Trigger/L1TMuonBayes/interface/Omtf/GoldenPattern.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFConfiguration.h>
#include <L1Trigger/L1TMuonBayes/interface/Omtf/OMTFinput.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "FWCore/MessageLogger/interface/MessageLogger.h"


PdfValueType GoldenPattern::pdfValue(unsigned int iLayer, unsigned int iRefLayer, unsigned int iBin, int refLayerPhiB) const {

  return pdfAllRef[iLayer][iRefLayer][iBin];
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
int GoldenPattern::meanDistPhiValue(unsigned int iLayer, unsigned int iRefLayer, int refLayerPhiB) const {
  //return meanDistPhi[iLayer][iRefLayer][0];
  return ( ( (meanDistPhi[iLayer][iRefLayer][1]*refLayerPhiB)>>myOmtfConfig->nPdfAddrBits() ) + meanDistPhi[iLayer][iRefLayer][0] );
  //assumes that the meanDistPhi[1] is float alpha from the fit to the phiB-phi distribution multiplied by 2^myOmtfConfig->nPdfAddrBits()
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
int GoldenPattern::propagateRefPhi(int phiRef, int etaRef, unsigned int iRefLayer){
  unsigned int iLayer = 2; //MB2 
  //if(etaRef>101) iLayer = 7;//RE2
  return phiRef + meanDistPhi[iLayer][iRefLayer][0];
  //FIXME if the meanDistPhiAlpha is non-zero, then meanDistPhi is alone not good for propagation of the phi
  //other value should be used, or the ref_layer phiB should be included
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPattern::normalisePdf(){

  unsigned int nPdfAddrBits = myOmtfConfig->nPdfAddrBits();
  unsigned int nPdfValBits = myOmtfConfig->nPdfValBits();
  double minP = myOmtfConfig->minPdfVal();
  for (unsigned int iRefLayer=0;iRefLayer<pdfAllRef[0].size();++iRefLayer){
    for (unsigned int iLayer=0;iLayer<pdfAllRef.size();++iLayer){
      double sum = 0.0;
      double value = 0.0;
      //Compute sum of probabilities
      for (unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){   
        value = pdfAllRef[iLayer][iRefLayer][iPdf];
	value -= (1<<nPdfValBits)-1;
	value *= log(minP)/((1<<nPdfValBits)-1);
	double tmp =  exp(-value);
	value = tmp - minP;
	sum+=value;  
      }      
      
      double pValLog = 0.0;
      int digitisedVal = 0;
      ///Normalise probabilites to 1.0, and convert to
      ///integer scale
      for (unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){
	double value = pdfAllRef[iLayer][iRefLayer][iPdf];
	value = pdfAllRef[iLayer][iRefLayer][iPdf];
	value -= (1<<nPdfValBits)-1;
	value *= log(minP)/(double)((1<<nPdfValBits)-1);
	value = exp(-value) - minP;
	value /= sum;
	pValLog = log(value);	
	digitisedVal = rint((1<<nPdfValBits)-1 - (pValLog/log(minP))*((1<<nPdfValBits)-1));
	pdfAllRef[iLayer][iRefLayer][iPdf] = digitisedVal;
      }
    }
  }  
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
std::ostream & operator << (std::ostream &out, const GoldenPattern & aPattern){

  out <<"GoldenPattern "<< aPattern.theKey <<std::endl;
  out <<"Number of reference layers: "<<aPattern.meanDistPhi[0].size()
          <<", number of measurement layers: "<<aPattern.pdfAllRef.size()
          <<std::endl;

  if(aPattern.meanDistPhi.empty()) return out;
  if(aPattern.pdfAllRef.empty()) return out;

  out<<"Mean dist phi per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer<<" (";
    for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
      for(unsigned int iPar = 0; iPar < aPattern.meanDistPhi[iLayer][iRefLayer].size(); iPar++)
        out<<std::setw(3)<<aPattern.meanDistPhi[iLayer][iRefLayer][iPar]<<"\t";
    }
    out<<")"<<std::endl;
  }

/*  if(aPattern.meanDistPhiCounts.size()){
    out<<"Counts number per layer:"<<std::endl;
    for (unsigned int iRefLayer=0;iRefLayer<aPattern.meanDistPhi[0].size();++iRefLayer){
      out<<"Ref layer: "<<iRefLayer<<" (";
      for (unsigned int iLayer=0;iLayer<aPattern.meanDistPhi.size();++iLayer){   
        out<<aPattern.meanDistPhiCounts[iLayer][iRefLayer]<<"\t";
      }
      out<<")"<<std::endl;
    }
  }*/

  unsigned int nPdfAddrBits = 7;
  out<<"PDF per layer:"<<std::endl;
  for (unsigned int iRefLayer=0;iRefLayer<aPattern.pdfAllRef[0].size();++iRefLayer){
    out<<"Ref layer: "<<iRefLayer;
    for (unsigned int iLayer=0;iLayer<aPattern.pdfAllRef.size();++iLayer){   
      out<<", measurement layer: "<<iLayer<<std::endl;
      for (unsigned int iPdf=0;iPdf<exp2(nPdfAddrBits);++iPdf){   
        out<<std::setw(2)<<aPattern.pdfAllRef[iLayer][iRefLayer][iPdf]<<" ";
      }
      out<<std::endl;
    }
  }

  return out;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
void GoldenPattern::reset() {
  for (unsigned int iLayer=0; iLayer < meanDistPhi.size(); ++iLayer) {
    for (unsigned int iRefLayer=0; iRefLayer < meanDistPhi[iLayer].size(); ++iRefLayer) {
      for (unsigned int iBin=0; iBin < meanDistPhi[iLayer][iRefLayer].size(); ++iBin) {
        meanDistPhi[iLayer][iRefLayer][iBin] = 0;
      }
    }
  }

  for (unsigned int iLayer=0; iLayer < distPhiBitShift.size(); ++iLayer) {
    for (unsigned int iRefLayer=0; iRefLayer < distPhiBitShift[iLayer].size(); ++iRefLayer) {
      distPhiBitShift[iLayer][iRefLayer] = 0;
    }
  }

  for (unsigned int iLayer=0; iLayer < pdfAllRef.size(); ++iLayer) {
    for (unsigned int iRefLayer=0; iRefLayer < pdfAllRef[iLayer].size(); ++iRefLayer) {
      for (unsigned int iPdf=0; iPdf < pdfAllRef[iLayer][iRefLayer].size(); ++iPdf) {
        pdfAllRef[iLayer][iRefLayer][iPdf] = 0;
      }
    }
  }

}
