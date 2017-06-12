#include <iostream>
#include <strstream>
#include <iomanip>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/plugins/OMTFHitAnalyzer.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfigMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "Math/VectorUtil.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include <TH2F.h>
#include "TFile.h"


OMTFHitAnalyzer::OMTFHitAnalyzer(const edm::ParameterSet& cfg):
theConfig(cfg),
g4SimTrackSrc(cfg.getParameter<edm::InputTag>("g4SimTrackSrc")){

  inputTokenDTPh = consumes<L1MuDTChambPhContainer>(theConfig.getParameter<edm::InputTag>("srcDTPh"));
  inputTokenDTTh = consumes<L1MuDTChambThContainer>(theConfig.getParameter<edm::InputTag>("srcDTTh"));
  inputTokenCSC = consumes<CSCCorrelatedLCTDigiCollection>(theConfig.getParameter<edm::InputTag>("srcCSC"));
  inputTokenRPC = consumes<RPCDigiCollection>(theConfig.getParameter<edm::InputTag>("srcRPC"));
  inputTokenSimHit = consumes<edm::SimTrackContainer>(theConfig.getParameter<edm::InputTag>("g4SimTrackSrc"));

  myInputMaker = new OMTFinputMaker();

  makeGoldenPatterns = theConfig.getParameter<bool>("makeGoldenPatterns");
  makeConnectionsMaps = theConfig.getParameter<bool>("makeConnectionsMaps");
  mergeXMLFiles = theConfig.getParameter<bool>("mergeXMLFiles");

  myOMTFConfig = 0;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFHitAnalyzer::~OMTFHitAnalyzer(){

  delete myOMTFConfig;
  //delete myOMTFConfigMaker;
  delete myOMTF;

}
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::configureProcesor(const OMTFConfiguration * omtfConfig,
    const L1TMuonOverlapParams* omtfPatterns, OMTFProcessor* omtfProc, unsigned int ptCode, int charge) {
  omtfProc->configure(myOMTFConfig);

  //myResults.assign(omtfConfig->nTestRefHits(),OMTFProcessor::resultsMap()); FIXME is it needed???

  const l1t::LUT* chargeLUT =  omtfPatterns->chargeLUT();
  const l1t::LUT* etaLUT =  omtfPatterns->etaLUT();
  const l1t::LUT* ptLUT =  omtfPatterns->ptLUT();

  unsigned int nGPs = omtfConfig->nGoldenPatterns();
  unsigned int address = 0;
  unsigned int iEta = 0;
  unsigned int iPt = ptCode;
  int iCharge = charge;
  /*for(unsigned int iGP=0; iGP < nGPs; ++iGP) {
    address = iGP;
    iEta = etaLUT->data(address);
    iCharge = chargeLUT->data(address)==0 ? -1:1;
    iPt = ptLUT->data(address);
    std::cout<<"iGP "<<iGP<<" etaLUT "<<iEta<<" iCharge "<<iCharge<<" iPt "<<iPt<<std::endl;
    if(iPt == ptCode && iCharge == charge) { //TODO check the ptCode
      Key aKey(iEta,iPt,iCharge,iGP);
      std::cout<<"adding GoldenPatternPdf4D "<<aKey<<std::endl;
      GoldenPatternPdf4D *aGP = new GoldenPatternPdf4D(aKey, omtfConfig);
      aGP->reset();
      omtfProc->addGP(aGP);
    }
  }*/

  Key aKey(iEta,iPt,iCharge, 0);
  std::cout<<"adding GoldenPatternPdf4D "<<aKey<<std::endl;
  GoldenPatternPdf4D *aGP = new GoldenPatternPdf4D(aKey, omtfConfig);
  aGP->reset();
  omtfProc->addGP(aGP);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

  const L1TMuonOverlapParamsRcd& omtfParamsRcd = iSetup.get<L1TMuonOverlapParamsRcd>();

  edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
  omtfParamsRcd.get(omtfParamsHandle);

  const L1TMuonOverlapParams* omtfParams = omtfParamsHandle.product();

  if (!omtfParams) {
    edm::LogError("L1TMuonOverlapTrackProducer") << "Could not retrieve parameters from Event Setup" << std::endl;
  }

  ///Initialise XML writer with default pdf.
  myWriter = new XMLConfigWriter(myOMTFConfig);

  ///For making the patterns use extended pdf width in phi, as pdf are later shifted by the mean value
  ///For low pt muons non shifted pdfs would go out of the default pdf range.
  L1TMuonOverlapParams omtfParamsMutable = *omtfParams;
  std::vector<int> generalParams = *omtfParamsMutable.generalParams();
  nPdfAddrBits = omtfParams->nPdfAddrBits();

  //if(!mergeXMLFiles) //the LUTs are then too big in the GoldenPatternPdf4D
  generalParams[L1TMuonOverlapParams::GENERAL_ADDRBITS] = nPdfAddrBits + 2; //2*nPdfAddrBits;
  omtfParamsMutable.setGeneralParams(generalParams);

  myOMTFConfig->configure(&omtfParamsMutable);
  int ptCode = theConfig.getParameter<int>("ptCode");
  int charge = theConfig.getParameter<int>("charge");
  configureProcesor(myOMTFConfig, omtfParams, myOMTF, ptCode, charge);
  //myOMTFConfigMaker = new OMTFConfigMaker(myOMTFConfig);

  std::cout<<"myOMTFConfig "<<*myOMTFConfig;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::beginJob(){

  myOMTFConfig = new OMTFConfiguration();
  myOMTF = new OMTFProcessor();
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
void OMTFHitAnalyzer::endJob(){
  std::cout<<"myOMTFConfig "<<*myOMTFConfig;

  if(makeGoldenPatterns && !makeConnectionsMaps) {
    myWriter->initialiseXMLDocument("OMTF");
    const std::map<Key, IGoldenPattern*> & myGPmap = myOMTF->getPatterns();
    for(auto itGP: myGPmap){
      if(!itGP.second->hasCounts()) continue;
      itGP.second->normalise(nPdfAddrBits);
    }


    ///Put back default value of the pdf width.
/*    L1TMuonOverlapParams omtfParamsMutable = *myOMTFConfig->getRawParams();
    std::vector<int> generalParams = *omtfParamsMutable.generalParams();
    generalParams[L1TMuonOverlapParams::GENERAL_ADDRBITS] = nPdfAddrBits;
    omtfParamsMutable.setGeneralParams(generalParams);
    myOMTFConfig->configure(&omtfParamsMutable);*/

    std::ostrstream fileName;
    fileName<<"bendinfDistr_ptCode"<<theConfig.getParameter<int>("ptCode")
        <<"_ch"<<theConfig.getParameter<int>("charge")<<".root";
    TFile* outfile = new TFile(fileName.str(), "RECREATE");
    for(auto itGP: myGPmap) {
      ////
      unsigned int iPt = theConfig.getParameter<int>("ptCode")+1;
      if(iPt>31)
        iPt = 200*2+1;
      else
        iPt = RPCConst::ptFromIpt(iPt)*2.0+1;//MicroGMT has 0.5 GeV step size, with lower bin edge  (uGMT_pt_code - 1)*step_size
      ////
      /*      if(itGP.first.thePtCode==iPt &&
          itGP.first.theCharge==theConfig.getParameter<int>("charge")) {
        //std::cout<<*itGP.second<<std::endl; FIXME
        myWriter->writeGPData(*((GoldenPattern*)(itGP.second)), dummyGP, dummyGP, dummyGP);
      }*/
      for(unsigned int iLayer = 0; iLayer<myOMTFConfig->nLayers(); ++iLayer) {
        for(unsigned int iRefLayer=0; iRefLayer<myOMTFConfig->nRefLayers(); ++iRefLayer) {
          std::ostrstream histName;

          histName<<"ipt_"<<itGP.first.thePtCode<<"_ch"<<itGP.first.theCharge<<"_layer_"<<iLayer<<"_refLayer_"<<iRefLayer<<" ";

          GoldenPatternPdf4D* gp4D = (static_cast<GoldenPatternPdf4D*>(itGP.second));
          unsigned int refLayerPhiBSize = gp4D->getPdf()[iLayer][iRefLayer].size();
          unsigned int layerPhiSize = (static_cast<GoldenPatternPdf4D*>(itGP.second))->getPdf()[iLayer][iRefLayer][0].size();
          cout<<"creating hist "<<histName.str()<<" refLayerPhiBSize "<<refLayerPhiBSize<<" layerPhiSize "<<layerPhiSize<<std::endl;

          if(refLayerPhiBSize == 1 ) {
            TH1F *h1 = new TH1F(histName.str(), histName.str(), layerPhiSize, -0.5, layerPhiSize-0.5);
            for(unsigned int iLayerPhi=0; iLayerPhi < layerPhiSize; iLayerPhi++) {
              h1->Fill(iLayerPhi, gp4D->getPdf()[iLayer][iRefLayer][0][iLayerPhi]);
            }
            h1->Write();
          }
          else {
            TH2F *h2 = new TH2F(histName.str(), histName.str(), refLayerPhiBSize, -0.5, refLayerPhiBSize-0.5,
                layerPhiSize, -0.5, layerPhiSize-0.5);

            for(unsigned int iRefLayerPhiB = 0; iRefLayerPhiB < refLayerPhiBSize; iRefLayerPhiB++) {
              for(unsigned int iLayerPhi=0; iLayerPhi < layerPhiSize; iLayerPhi++) {
                h2->Fill(iRefLayerPhiB, iLayerPhi, gp4D->getPdf()[iLayer][iRefLayer][iRefLayerPhiB][iLayerPhi]);
              }
            }
            h2->Write();
          }
        }
      }
    }
    outfile->Close();
  }


}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::writeMergedGPs(){

/*  const std::map<Key, IGoldenPattern*> & myGPmap = myOMTF->getPatterns();

  GoldenPattern *dummy = new GoldenPattern(Key(0,0,0), myOMTFConfig);
  dummy->reset();

  unsigned int iPtMin = 9;
  Key aKey = Key(0, iPtMin, 1);
  while(myGPmap.find(aKey)!=myGPmap.end()){

    GoldenPattern *aGP1 = myGPmap.find(aKey)->second;
    GoldenPattern *aGP2 = dummy;
    GoldenPattern *aGP3 = dummy;
    GoldenPattern *aGP4 = dummy;

    ++aKey.thePtCode;
    while(myGPmap.find(aKey)==myGPmap.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
    if(aKey.thePtCode<=401 && myGPmap.find(aKey)!=myGPmap.end()) aGP2 =  myGPmap.find(aKey)->second;

    if(aKey.thePtCode>71){
      ++aKey.thePtCode;
      while(myGPmap.find(aKey)==myGPmap.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
      if(aKey.thePtCode<=401 && myGPmap.find(aKey)!=myGPmap.end()) aGP3 =  myGPmap.find(aKey)->second;

      ++aKey.thePtCode;
      while(myGPmap.find(aKey)==myGPmap.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
      if(aKey.thePtCode<=401 && myGPmap.find(aKey)!=myGPmap.end()) aGP4 =  myGPmap.find(aKey)->second;
    }
    ++aKey.thePtCode;
    while(myGPmap.find(aKey)==myGPmap.end() && aKey.thePtCode<=401) ++aKey.thePtCode;    
    myWriter->writeGPData(*aGP1,*aGP2, *aGP3, *aGP4);

    ///Write the opposite charge.
    Key aTmpKey = aGP1->key();
    aTmpKey.theCharge = -1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP1 =  myGPmap.find(aTmpKey)->second;
    else aGP1 = dummy;

    aTmpKey = aGP2->key();
    aTmpKey.theCharge = -1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP2 =  myGPmap.find(aTmpKey)->second;
    else aGP2 = dummy;

    aTmpKey = aGP3->key();
    aTmpKey.theCharge = -1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP3 =  myGPmap.find(aTmpKey)->second;
    else aGP3 = dummy;

    aTmpKey = aGP4->key();
    aTmpKey.theCharge = -1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP4 =  myGPmap.find(aTmpKey)->second;
    else aGP4 = dummy;

    myWriter->writeGPData(*aGP1,*aGP2, *aGP3, *aGP4);
  }*/
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup){

  if(mergeXMLFiles) return;

  ///Get the simulated muon parameters
  const SimTrack* aSimMuon = findSimMuon(iEvent,evSetup);
  if(!aSimMuon){
    edm::LogError("OMTFHitAnalyzer")<<"No SimMuon found in the event!";
    return;
  }

  //std::cout<<"new event ";
  //cout<<"aSimMuon->momentum().pt "<<aSimMuon->momentum().pt()<<std::endl;

  myInputMaker->initialize(evSetup, myOMTFConfig);

  edm::Handle<L1MuDTChambPhContainer> dtPhDigis;
  edm::Handle<L1MuDTChambThContainer> dtThDigis;
  edm::Handle<CSCCorrelatedLCTDigiCollection> cscDigis;
  edm::Handle<RPCDigiCollection> rpcDigis;

  ///Filter digis by dropping digis from selected (by cfg.py) subsystems
  if(!theConfig.getParameter<bool>("dropDTPrimitives")){
    iEvent.getByToken(inputTokenDTPh,dtPhDigis);
    iEvent.getByToken(inputTokenDTTh,dtThDigis);
  }
  if(!theConfig.getParameter<bool>("dropRPCPrimitives"))
    iEvent.getByToken(inputTokenRPC,rpcDigis);
  if(!theConfig.getParameter<bool>("dropCSCPrimitives"))
    iEvent.getByToken(inputTokenCSC,cscDigis);

  //l1t::tftype mtfType = l1t::tftype::bmtf;
  l1t::tftype mtfType = l1t::tftype::omtf_pos;
  //l1t::tftype mtfType = l1t::tftype::emtf_pos;

  ///Loop over all processors, each covering 60 deg in phi
  for(unsigned int iProcessor=0;iProcessor<6;++iProcessor) {

    ///Input data with phi ranges shifted for each processor, so it fits 11 bits range
    OMTFinput myInput = myInputMaker->buildInputForProcessor(dtPhDigis.product(),
        dtThDigis.product(),
        cscDigis.product(),
        rpcDigis.product(),
        iProcessor,
        mtfType);


    if(makeGoldenPatterns)
      myOMTF->fillCounts(iProcessor, myInput, aSimMuon);


    std::ostrstream ostr;
    bool wasHit = false;
    for(unsigned int iLayer=0; iLayer < myOMTFConfig->nLayers(); ++iLayer){
      const OMTFinput::vector1D& layerHits = myInput.getLayerData(iLayer, false);
      ostr<<"layer "<<iLayer<<" ";
      for(unsigned int i = 0; i < layerHits.size(); i++) {
        if(layerHits[i] != 5400) {
          wasHit = true;
          ostr<<std::setw(4)<<layerHits[i]<<" ";
        }
        else
          ostr<<"     ";
      }
      ostr<<"\n";
    }
    if(wasHit) {
      //std::cout<<"iProcessor "<<iProcessor<<std::endl<<ostr.str()<<std::endl;
    }

  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
const SimTrack * OMTFHitAnalyzer::findSimMuon(const edm::Event &ev, const edm::EventSetup &es, const SimTrack * previous){

  const SimTrack * result = 0;
  edm::Handle<edm::SimTrackContainer> simTks;
  ev.getByToken(inputTokenSimHit,simTks);

  for (std::vector<SimTrack>::const_iterator it=simTks->begin(); it< simTks->end(); it++) {
    const SimTrack & aTrack = *it;
    if ( !(aTrack.type() == 13 || aTrack.type() == -13) )continue;
    if(previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(),previous->momentum())<0.07) continue;
    if ( !result || aTrack.momentum().pt() > result->momentum().pt()) result = &aTrack;
  }
  return result;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(OMTFHitAnalyzer);
