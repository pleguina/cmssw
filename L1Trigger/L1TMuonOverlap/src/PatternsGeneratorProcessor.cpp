/*
 * PatternsGeneratorProcessor.cpp
 *
 *  Created on: Aug 3, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PatternsGeneratorProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternParametrised.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include <TLinearFitter.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include <climits>

PatternsGeneratorProcessor::~PatternsGeneratorProcessor() {
  // TODO Auto-generated destructor stub
}

void PatternsGeneratorProcessor::generatePatterns() {
	path = "/afs/cern.ch/work/k/kpijanow/public/CMSSW_9_2_0/src/L1Trigger/L1TMuonOverlap/test/expert/";
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
	std::cout<<"#####   gneratePatterns   #####"<<std::endl;

	for(auto& it : theGPs){
		it->reset();
	}

	//TH1::AddDirectory(kFALSE);
	getDTLayers();

	for(unsigned int iGroup = 0; iGroup < mergedPartters.size(); iGroup++)
	{
		///Mean dist phi data
		for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
		  for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer)
		  {
          //check if refLayer is a DT
          if(std::find(refLayersDT.begin(), refLayersDT.end(), iRefLayer) != refLayersDT.end()){
            //histograms of DT layers vs DT_bending
            if(myOmtfConfig->getRefToLogicNumber().at(iRefLayer) + 1 == (int)iLayer){
              generateDTSelfPattern(iGroup, iRefLayer, iLayer);
            }
            else{
              //generateDTPattern(iGroup, iRefLayer, iLayer);
              copyDTPattern(iGroup, iRefLayer, iLayer);
            }
          }
          else{
            generateNonDTPattern(iGroup, iRefLayer, iLayer);
          }
			  //}
		  }
		}
		std::cout<<"######\t\tGROUP"<<iGroup<<"\t\t######"<<std::endl;
	}
	std::cout<<"######\t\t END \t\t######"<<std::endl;
}

void PatternsGeneratorProcessor::generateDTPattern(int iGroup, int iRefLayer, int iLayer){
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
	int iMerge = mergedPartters[iGroup].size();
	std::vector<TH2F*> h;
	bool empty = true;

	for(int k = 0; k < iMerge; k++)
	{
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

		h.push_back(getHistogram2D(ipt, charge,
				iLayer, iRefLayer));
		TH2F* hRef = getHistogram2D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
    if(h.at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is empty
      hRef->Delete();
      continue;
    }

		double sigY = getMeanSigma(h.at(k)) * 0.8; 		//sigma in Y direction, multiplication by 0.8 is empirical
		double norm =  h.at(k)->GetSum()/hRef->GetSum(); //1;//

		const float logC = log(norm * 1.0 / sqrt(2 * TMath::Pi()) / sigY);
		const float minPlog =  log(myOmtfConfig->minPdfVal());
		const int nPdfValBits = myOmtfConfig->nPdfValBits();

		///Pdf data
		for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
		{
			float pVal = logC - (iPdf - (int)(1<<myOmtfConfig->nPdfAddrBits())/2) * (iPdf - (int)(1<<myOmtfConfig->nPdfAddrBits())/2)/2.0/sigY/sigY;
			if(pVal<minPlog || hRef->GetSum() == 0){
				theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
			  continue;
			}

			int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
			theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
		}
		hRef->Delete();
		empty = false;
	}

	if(!empty)
	{
    //merging histograms
	  TH2F* hMerge = new TH2F();
	  unsigned int i  = 0;
	  //find non empty histogram
	  while (h.size() > i && !hMerge->GetSum())
	  {
	    if(h.at(i)->GetSum()) {
	      hMerge->Delete();
	      hMerge = (TH2F*)h.at(0)->Clone("hMerge");
	    }
	    i++;
	  }
    hMerge->Reset();
    for(int i = 0; i < iMerge; i++)
    {
      if(h.at(i)->GetSum() && h.at(i)->GetXaxis()->GetXmax() == hMerge->GetXaxis()->GetXmax())
        hMerge->Add(h.at(i));
    }

    //fitting
    TF1 *f1 = new TF1("f1","pol1");
    hMerge->Fit("f1", "QN");

    //function: y = A * x + B
    double B = f1->GetParameter(0);
    double A = f1->GetParameter(1);

    //saving calculated values to theGPs
    for(int k = 0; k < iMerge; k++)
    {
      int value = rint(B) - h.at(k)->GetNbinsX()/2 + h.at(k)->GetNbinsX()/2 * A;
      theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(value, iLayer, iRefLayer, 0);
      value = rint(A * pow(2, myOmtfConfig->nPdfAddrBits())); //multiply A by power of 2 to save it as integer
      theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(value, iLayer, iRefLayer, 1);
    }
    hMerge->Delete();
    f1->Delete();
	}
	for(auto &a : h)
	  a->Delete();
	h.clear();
}

void PatternsGeneratorProcessor::generateNonDTPattern(int iGroup, int iRefLayer, int iLayer){
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
	int iMerge = mergedPartters[iGroup].size();
	std::vector<TH1F*> h;
	int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);

  //save PDF's
	for(int k = 0; k < iMerge; k++)
	{
		unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
		int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

		TH1F* hRef = getHistogram1D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
		if(h.at(k)->GetSum() < emptyHistogramCount) //if this happens then histogram is empty
    {
		  hRef->Delete();
      continue;
    }

		theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - h.at(k)->GetNbinsX()/2, iLayer, iRefLayer, 0);

		const float minPlog =  log(myOmtfConfig->minPdfVal());
		const int nPdfValBits = myOmtfConfig->nPdfValBits();
		///Pdf data
		for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
		{
			// +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
			float pVal = log(h.at(k)->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/hRef->GetSum()); //h.at(k)->GetSum()); //
			if(pVal < minPlog || hRef->GetSum() == 0){
				theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
			  continue;
			}
			int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
			theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
		}
		hRef->Delete();
	}

  for(auto &a : h)
    a->Delete();
  h.clear();
}

void PatternsGeneratorProcessor::getDTLayers(){
	const std::map<int, int> hwToLogicLayer = myOmtfConfig->getHwToLogicLayer();
	layersDT.clear();
	refLayersDT.clear();

	layersDT.push_back(hwToLogicLayer.at(101));
	layersDT.push_back(hwToLogicLayer.at(102));
	layersDT.push_back(hwToLogicLayer.at(103));

	for(auto& L : layersDT){
		for(unsigned int i = 0; i < myOmtfConfig->getRefToLogicNumber().size(); i++){
			if(myOmtfConfig->getRefToLogicNumber().at(i) == L)
				refLayersDT.push_back(i);
		}
	}
}

void PatternsGeneratorProcessor::generateDTSelfPattern(int iGroup, int iRefLayer, int iLayer){
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
  int iMerge = mergedPartters[iGroup].size();
  std::vector<TH2F*> h;
  int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);

  for(int k = 0; k < iMerge; k++)
  {
    TH1D* hProj = (TH1D*)h.at(k)->ProjectionY("hProj");

    theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - hProj->GetNbinsX()/2, iLayer, iRefLayer, 0);

    const float minPlog =  log(myOmtfConfig->minPdfVal());
    const int nPdfValBits = myOmtfConfig->nPdfValBits();

    ///Pdf data
    for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
    {
      // +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
      float pVal = log(hProj->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/h.at(k)->GetSum());
      if(pVal < minPlog || hProj->GetSum() == 0){
        theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
        continue;
      }

      int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
      theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
    }

    delete hProj;
  }

  for(auto &a : h)
    delete a;
  h.clear();
}

void PatternsGeneratorProcessor::copyDTPattern(int iGroup, int iRefLayer, int iLayer){
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
  int iMerge = mergedPartters[iGroup].size();
  int histCount = 0;
  std::vector<TH2F*> h;
  int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);
  meanDistPhi = 0;
  for(auto a: h){
    if(a->GetSum() > emptyHistogramCount){
      meanDistPhi += a->ProjectionY("hProj")->GetMean();
      histCount++;}
  }
  if(histCount)
    meanDistPhi = rint(meanDistPhi / histCount);

  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

    TH2F* hRef = getHistogram2D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
    if(h.at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is empty
      delete hRef;
      continue;
    }
    TH1D* hProj = (TH1D*)h.at(k)->ProjectionY("hProj");
    theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - h.at(k)->GetNbinsX()/2, iLayer, iRefLayer, 0);
    const float minPlog =  log(myOmtfConfig->minPdfVal());
    const int nPdfValBits = myOmtfConfig->nPdfValBits();
    ///Pdf data
    for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
    {
      // +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
      float pVal = log(hProj->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/hRef->GetSum()); //h.at(k)->GetSum()); //
      if(pVal < minPlog || hRef->GetSum() == 0){
        theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
        continue;
      }
      int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
      if(ipt == 11 && iRefLayer == 5 && iLayer == 0 && charge == 1)
        std::cout<<iPdf<<" "<<digitisedVal<<" "<<minPlog<<" "<<pVal<<" "<<hRef->GetSum()<<" "<<hProj->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)<<std::endl;
      theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
    }
    delete hRef;
    delete hProj;
  }

  for(auto &a : h)
    delete a;
  h.clear();
}

int PatternsGeneratorProcessor::getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH2F*>* h)
{
  int meanDistPhi = 0;
  int meanDistCount = 0;
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
  int iMerge = mergedPartters[iGroup].size();
  //read histograms form a group and calculate meanDistPhi
  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;
    (*h).push_back(getHistogram2D(ipt, charge, iLayer, iRefLayer));
    if((*h).at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is considered empty
      continue;
    }
    if((ipt == 9 || ipt == 10) && iRefLayer == 5 && iLayer == 0 && charge == 1)
      std::cout<<"MEANDISTPHI "<<(*h).at(k)->GetMean()<<" "<<(*h).at(k)->GetSum()<<" "<<ipt<<std::endl;
    meanDistPhi += (*h).at(k)->GetMean();//* (*h).at(k)->GetSum();
    meanDistCount++;//= (*h).at(k)->GetSum();
  }
//  if(meanDistCount)
//    meanDistPhi = rint(meanDistPhi/meanDistCount);
  if(meanDistCount)
    meanDistPhi = meanDistPhi / meanDistCount;

  return meanDistPhi;
}

int PatternsGeneratorProcessor::getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH1F*>* h)
{
  int meanDistPhi = 0;
  int meanDistCount = 0;
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPatterns();
  int iMerge = mergedPartters[iGroup].size();
  //read histograms form a group and calculate meanDistPhi
  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;
    (*h).push_back(getHistogram1D(ipt, charge, iLayer, iRefLayer));
    if((*h).at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is considered empty
      continue;
    }
    meanDistPhi += (*h).at(k)->GetMean();//* (*h).at(k)->GetSum();
    meanDistCount++;//= (*h).at(k)->GetSum();
  }
//  if(meanDistCount)
//    meanDistPhi = rint(meanDistPhi/meanDistCount);
  if(meanDistCount)
    meanDistPhi = meanDistPhi / meanDistCount;

  return meanDistPhi;
}

TH1F* PatternsGeneratorProcessor::getHistogram1D(int ipt, int charge, int iLayer, int iRefLayer)
{
  std::ostringstream s1;
  s1 << path << "bendingDistr_ptCode_" << ipt << "_ch_" << charge << ".root";
  TH1F* h = new TH1F();

  if(access( (s1.str()).c_str(), F_OK ) == -1){
    std::cout<<"NO ACCESS to "<<s1.str().c_str()<<std::endl;
    return h;
  }

  TFile* file = TFile::Open((TString)s1.str());
  h->Delete();
  h = (TH1F*)file->Get(Form("ipt_%i_ch%i_layer_%i_refLayer_%i", ipt, charge,
      iLayer, iRefLayer));
  h->SetDirectory(0);
  delete file;
  return h;
}

TH2F* PatternsGeneratorProcessor::getHistogram2D(int ipt, int charge, int iLayer, int iRefLayer)
{
  std::ostringstream s1;
  s1 << path << "bendingDistr_ptCode_" << ipt << "_ch_" << charge << ".root";
  TH2F* h = new TH2F();

  if(access( (s1.str()).c_str(), F_OK ) == -1){
    std::cout<<"NO ACCESS to "<<s1.str().c_str()<<std::endl;
    return h;
  }

  TFile* file = TFile::Open((TString)s1.str());
  h->Delete();
  h = (TH2F*)file->Get(Form("ipt_%i_ch%i_layer_%i_refLayer_%i", ipt, charge,
      iLayer, iRefLayer));
  h->SetDirectory(0);
  delete file;
  return h;
}

double PatternsGeneratorProcessor::getMeanSigma(TH2F* h)
{
  //mean sigma from slices
  TObjArray aSlices;
  h->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);

  double meanSigma = 0;
  double denom = 0;
  TH1D * as = (TH1D*)aSlices[2]->Clone("as");
  for(int d = 0; d < as->GetNbinsX(); d++)
  {
    //error of the value of sigma cannot be grater than 20%
    if(as->GetBinError(d)/as->GetBinContent(d) < 0.2 && as->GetBinError(d) > 0 && as->GetBinContent(d) > 0)
    {
      meanSigma += as->GetBinContent(d) * 1/as->GetBinError(d);
      denom += 1/as->GetBinError(d);
    }
    else
    {
      as->SetBinContent(d, 0);
      as->SetBinError(d, 0);
    }
  }

  if(denom > 0)
    meanSigma = meanSigma / denom;

  //if sigma is lower than 1 bin than all values are in one bin ie. it is a histogram of DT vs itself
  //the value was chosen to match PDF values from Arthur's patterns
  if (meanSigma < 0.25)
    meanSigma = 0.25;

  delete as;
  return meanSigma;
}

void PatternsGeneratorProcessor::generateThresholds() {
/*  for(auto& gp : theGPs) {
    std::vector<unsigned int> thresholds(gp->getPdf()[0].size(), 0);
    for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      std::vector<int> maxPdfValues(gp->getPdf()[0].size());
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
        maxPdfValues[iRefLayer] = *std::max_element(gp->getPdf()[iLayer][iRefLayer].begin(), gp->getPdf()[iLayer][iRefLayer].end() );
        int pdfVal = 0;//maxPdfValues[iRefLayer]
        gp->setPdfValue(0, iLayer, iRefLayer, 0);//for the version with thresholds in the bin 0 of the pdf we  keep the max vale
        //it is used when given layer is not fired
      }

      int  smallestMax  = INT_MAX;
      unsigned int smallestMaxPos = maxPdfValues.size();
      for(unsigned int i = 0; i < maxPdfValues.size(); i++) {
        if(maxPdfValues[i] > 0 && smallestMax > maxPdfValues[i]) {
          smallestMax = maxPdfValues[i];
          smallestMaxPos = i;
        }
      }
      if(smallestMax == INT_MAX) {
        smallestMax = 0;
        smallestMaxPos = maxPdfValues.size();
      }

      thresholds[iRefLayer] = smallestMax;
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
        if(iLayer == smallestMaxPos)
          continue;

        int minPdfVal = INT_MAX;//todo
        for(auto& pdfVal : gp->getPdf()[iLayer][iRefLayer]) {
          if( (pdfVal > 0) && (minPdfVal > pdfVal) )
            minPdfVal = pdfVal;
        }
        if(minPdfVal != INT_MAX)
          thresholds[iRefLayer] += minPdfVal;
        //std::cout<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" minPdfVal "<<minPdfVal<<" thr "<<thresholds[iRefLayer]<<std::endl;
      }
      if(gp->key().thePt <= 10) {
        thresholds[iRefLayer] = 0;
      }
      else if(gp->key().thePt >= 50) {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 1.35;
      }
      else if(gp->key().thePt >= 33) {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 1.35;
      }
      else {
        thresholds[iRefLayer] =  thresholds[iRefLayer] * 0.4;
      }
    }
    gp->setThresholds(thresholds);
  }*/
}


void PatternsGeneratorProcessor::modifyPdfs() {
  for(auto& gp : theGPs) {
    GoldenPatternParametrised* newGp = new GoldenPatternParametrised(gp.get());
    delete newGp;
    break;

    /*if(gp->key().thePt <= 25 || gp->key().thePt >= 100)
      continue;
    for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
        if(iRefLayer == 1 || iRefLayer == 3 || iRefLayer == 4) {
          if(iLayer == 6 || iLayer == 7 || iLayer == 8 || iLayer == 15 || iLayer == 16 || iLayer == 17) {
            for(unsigned int iPdf = 1; iPdf < gp->getPdf().at(iLayer).at(iRefLayer).size(); iPdf++) {
              int newPdfVal = (gp->pdfValue(iLayer, iRefLayer, iPdf) - 25) * 1.5;
              if(newPdfVal < 0)
                newPdfVal =  0;
              gp->setPdfValue(newPdfVal, iLayer, iRefLayer, iPdf);

            }
            for(unsigned int iPdf = 1; iPdf < gp->getPdf().at(iLayer).at(iRefLayer).size() -1; iPdf++) {
              if(gp->pdfValue(iLayer, iRefLayer, iPdf) != 0 && gp->pdfValue(iLayer, iRefLayer, iPdf) <= 30 &&
                 gp->pdfValue(iLayer, iRefLayer, iPdf+1) == 0 ) {
                gp->setPdfValue(0, iLayer, iRefLayer, iPdf);
                break;
              }
            }

            for(unsigned int iPdf = gp->getPdf().at(iLayer).at(iRefLayer).size() -1; iPdf >= 1 ; iPdf--) {
              if(gp->pdfValue(iLayer, iRefLayer, iPdf) != 0 && gp->pdfValue(iLayer, iRefLayer, iPdf) <= 30 &&
                 gp->pdfValue(iLayer, iRefLayer, iPdf-1) == 0 ) {
                gp->setPdfValue(0, iLayer, iRefLayer, iPdf);
                break;
              }
            }
          }
        }
      }
    }*/
  }
}

void PatternsGeneratorProcessor::modifyPdfs1() {
  /*for(auto& gp : theGPs) {
    for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[0].size(); ++iRefLayer) {
      for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {

        if(gp->getPdf()[iLayer][iRefLayer][0] > 0 || gp->getPdf()[iLayer][iRefLayer][1] > 0 ||
           gp->getPdf()[iLayer][iRefLayer][127] > 0 || gp->getPdf()[iLayer][iRefLayer][126] > 0 ) {
          std::cout<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<std::endl;
          for(unsigned int iPdf = 1; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
            std::cout<<iPdf<<"\t"<<gp->getPdf()[iLayer][iRefLayer][iPdf]<<std::endl;
          }
          std::cout<<std::endl;
        }
      }
    }
  }*/

  for(auto& gp : theGPs) {
    if(gp->key().thePt >= 14 || gp->key().thePt == 0)
      continue;
    for(unsigned int iLayer = 0; iLayer < gp->getPdf().size(); ++iLayer) {
      for(unsigned int iRefLayer=0; iRefLayer < gp->getPdf()[iLayer].size(); ++iRefLayer) {
        //if(iRefLayer == 1 || iRefLayer == 3 || iRefLayer == 4)
        {
          if( ( (gp->key().thePt == 9 || gp->key().thePt == 10) && (iLayer == 1 || iLayer == 3 || (iLayer == 5 && iRefLayer == 5) ) ) ||
              ( (gp->key().thePt == 11|| gp->key().thePt == 13) && (iLayer == 1 ) )    ) {
            std::cout<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<std::endl;
            gp->setDistPhiBitShift(1, iLayer, iRefLayer);
            boost::multi_array<omtfPdfValueType, 1> pdfBuf(boost::extents[gp->getPdf()[iLayer][iRefLayer].size()]);
            pdfBuf = gp->getPdf()[iLayer][iRefLayer];
            for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
              assert(gp->getPdf()[iLayer][iRefLayer][iPdf] == pdfBuf[iPdf]);
              gp->setPdfValue(0, iLayer, iRefLayer, iPdf);
            }
            for(unsigned int iBinNew = 32; iBinNew < pdfBuf.size() - 32; iBinNew++) {
              int iBinOld = (iBinNew - 32) * 2;
              double pdfVal = (pdfBuf[iBinOld] + pdfBuf[iBinOld + 1]) / 2.;
              gp->setPdfValue(rint(pdfVal), iLayer, iRefLayer, iBinNew);
            }

            TLinearFitter linearFitter;
            linearFitter.StoreData(kTRUE);
            std::ostringstream ostr;
            ostr<<"fit_lay_"<<iLayer<<"_refLay_"<<iRefLayer;
            TF1 linearFit("linearFit", "pol2", 0, gp->getPdf()[iLayer][iRefLayer].size());
            linearFitter.SetFormula(linearFit.GetFormula());

            for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
              if(gp->pdfValue(iLayer, iRefLayer, iPdf) > 0) {
                double x = iPdf;
                linearFitter.AddPoint(&x, gp->pdfValue(iLayer, iRefLayer, iPdf));
              }
            }
            int res = linearFitter.Eval();
            std::cout<<gp->key()<<" iLayer "<<iLayer<<" iRefLayer "<<iRefLayer<<" linearFitter result "<<res<<std::endl;
            if(res == 0) {
              for(unsigned int iPdf = 0; iPdf < gp->getPdf()[iLayer][iRefLayer].size(); iPdf++) {
                if(gp->getPdf()[iLayer][iRefLayer][iPdf] != 0)
                  continue;
                double pdfVal = 0;
                unsigned int iBinPow = 1;
                for(unsigned int i = 0; i < 3; i++) {
                  //std::cout<<linearFit.GetParName(i)<<" "<<linearFit.GetParameter(i)<<" ";

                  pdfVal += linearFit.GetParameter(i) * iBinPow ;
                  iBinPow *= iPdf;
                }
                if(pdfVal > 0) {
                  gp->setPdfValue(rint(pdfVal), iLayer, iRefLayer, iPdf);
                }
              }
            }
          }
        }
      }
    }
  }
}

