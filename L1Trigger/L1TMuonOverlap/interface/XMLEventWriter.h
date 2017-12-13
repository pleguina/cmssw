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

  virtual void observeProcesorEmulation(unsigned int iProcessor, l1t::tftype mtfType,  const OMTFinput &input,
      const std::vector<AlgoMuon>& algoCandidates,
      std::vector<AlgoMuon>& gbCandidates,
      const std::vector<l1t::RegionalMuonCand> & candMuons);

  virtual void observeEventBegin(const edm::Event& iEvent);

  virtual void endJob();
private:
  const OMTFConfiguration* omtfConfig;
  XMLConfigWriter xmlWriter;
  xercesc::DOMElement* topElement;

  std::string fName;

  unsigned int eventNum = 0;
};

#endif /* INTERFACE_XMLEVENTWRITER_H_ */
