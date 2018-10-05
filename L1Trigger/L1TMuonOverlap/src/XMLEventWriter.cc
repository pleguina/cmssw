/*
 * XMLEventWriter.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/XMLEventWriter.h"

XMLEventWriter::XMLEventWriter(const OMTFConfiguration* aOMTFConfig, std::string fName):
omtfConfig(aOMTFConfig), xmlWriter(aOMTFConfig), currentElement(nullptr), fName(fName) {
  //std::string fName = "OMTF";
  xmlWriter.initialiseXMLDocument("OMTF");
  eventNum = 0;
};


XMLEventWriter::~XMLEventWriter() {
  // TODO Auto-generated destructor stub
}

void XMLEventWriter::observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
    const std::vector<AlgoMuon>& algoCandidates,
    std::vector<AlgoMuon>& gbCandidates,
    const std::vector<l1t::RegionalMuonCand> & candMuons)
{
  if(eventNum > 1000)
    return;

  int endcap =  (mtfType == l1t::omtf_neg) ? -1 : ( ( mtfType == l1t::omtf_pos) ? +1 : 0 );
  OmtfName board(iProcessor, endcap);

  if(candMuons.empty())
    return;

  //if(currentElement == nullptr)
    currentElement = xmlWriter.writeEventHeader(eventId);

  xercesc::DOMElement * aProcElement = xmlWriter.writeEventData(currentElement, board, input);
  for(unsigned int iRefHit=0;iRefHit < omtfConfig->nTestRefHits();++iRefHit){
    ///Dump only regions, where a candidate was found
    const AlgoMuon& algoMuon = algoCandidates.at(iRefHit);//charge=0 means ignore charge
    if(algoMuon.isValid()) {
      xmlWriter.writeAlgoMuon(aProcElement,iRefHit,algoMuon);
      /*if(dumpDetailedResultToXML){
        for(auto & itKey: results[iRefHit])
          xmlWriter.writeResultsData(aProcElement, iRefHit, itKey.first,itKey.second);
      }*/
    }
  }
  for (auto & candMuon :  candMuons) xmlWriter.writeCandMuon(aProcElement, candMuon);

}

void XMLEventWriter::observeEventBegin(const edm::Event& iEvent) {
  eventNum++;
  if(eventNum > 1000) //due to some bug if more events is written the memory consumption s very big and program crashes
    return;
  //currentElement = xmlWriter.writeEventHeader(iEvent.id().event());
  eventId = iEvent.id().event();
}

void XMLEventWriter::observeEventEnd(const edm::Event& iEvent) {
  currentElement = nullptr;
}

void XMLEventWriter::endJob() {
  xmlWriter.finaliseXMLDocument(fName);
}

