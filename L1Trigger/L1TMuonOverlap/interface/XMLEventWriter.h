/*
 * XMLEventWriter.h
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#ifndef INTERFACE_XMLEVENTWRITER_H_
#define INTERFACE_XMLEVENTWRITER_H_

#include "L1Trigger/L1TMuonOverlap/interface/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"

class XMLEventWriter: public IOMTFEmulationObserver {
public:
  XMLEventWriter(const OMTFConfiguration* aOMTFConfig, std::string fName);

  virtual ~XMLEventWriter();

  void observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
      const AlgoMuons& algoCandidates,
      const AlgoMuons& gbCandidates,
      const std::vector<l1t::RegionalMuonCand> & candMuons) override;

  void observeEventBegin(const edm::Event& iEvent) override;

  void observeEventEnd(const edm::Event& iEvent) override;

  void endJob() override;
private:
  const OMTFConfiguration* omtfConfig;
  XMLConfigWriter xmlWriter;
  xercesc::DOMElement* currentElement;

  std::string fName;

  unsigned int eventNum = 0;

  unsigned int eventId = 0;
};

#endif /* INTERFACE_XMLEVENTWRITER_H_ */
