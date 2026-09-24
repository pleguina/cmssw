/*
 * XMLEventWriter.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Common/interface/EventBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfName.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/XMLEventWriter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <boost/property_tree/xml_parser.hpp>

#include <bitset>

XMLEventWriter::XMLEventWriter(const OMTFConfiguration* aOMTFConfig, std::string fName, int eventsPerFile)
    : omtfConfig(aOMTFConfig), fName(fName), eventsPerFile(eventsPerFile) {
  eventNum = 0;
  initTree();
};

void XMLEventWriter::initTree() {
  tree.clear();
  unsigned int version = omtfConfig->patternsVersion();
  const unsigned int mask16bits = 0xFFFF;
  version &= mask16bits;
  std::ostringstream stringStr;
  stringStr << "0x" << std::hex << std::setfill('0') << std::setw(4) << version;
  tree.put("OMTF.<xmlattr>.version", stringStr.str());
}

void XMLEventWriter::flushCurrentTreeToFile() {
  // Build part filename: <base>_part000<ext>
  std::string baseName = fName;
  std::string ext;
  auto dotPos = fName.rfind('.');
  if (dotPos != std::string::npos) {
    baseName = fName.substr(0, dotPos);
    ext = fName.substr(dotPos);
  }
  std::ostringstream ss;
  ss << baseName << "_part" << std::setfill('0') << std::setw(3) << fileIndex << ext;
  std::string partFileName = ss.str();

  edm::LogInfo("l1tOmtfEventPrint") << "XMLEventWriter: writing part file " << partFileName;
  boost::property_tree::write_xml(
      partFileName, tree, std::locale(), boost::property_tree::xml_parser::xml_writer_make_settings<std::string>(' ', 2));

  fileIndex++;
  initTree();  // reset tree for next batch
}

XMLEventWriter::~XMLEventWriter() {}

void XMLEventWriter::observeProcesorBegin(unsigned int iProcessor, l1t::tftype mtfType) {
  procTree.clear();

  int endcap = (mtfType == l1t::omtf_neg) ? -1 : ((mtfType == l1t::omtf_pos) ? +1 : 0);
  OmtfName board(iProcessor, endcap, omtfConfig);
  procTree.add("<xmlattr>.board", board.name());
  procTree.add("<xmlattr>.iProcessor", iProcessor);

  std::ostringstream stringStr;
  stringStr << (board.position() == 1 ? "+" : "") << board.position();
  procTree.add("<xmlattr>.position", stringStr.str());
}

void XMLEventWriter::observeProcesorEmulation(unsigned int iProcessor,
                                              l1t::tftype mtfType,
                                              const std::shared_ptr<OMTFinput>& input,
                                              const AlgoMuons& algoCandidates,
                                              const AlgoMuons& gbCandidates,
                                              const FinalMuons& finalMuons) {
  int endcap = (mtfType == l1t::omtf_neg) ? -1 : ((mtfType == l1t::omtf_pos) ? +1 : 0);
  OmtfName board(iProcessor, endcap, omtfConfig);

  if (finalMuons.empty())
    return;

  // Layer/Hit generation removed per user request
  /*
  for (unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) {
    boost::property_tree::ptree layerTree;

    for (unsigned int iHit = 0; iHit < input->getMuonStubs()[iLayer].size(); ++iHit) {
      int hitPhi = input->getPhiHw(iLayer, iHit);
      if (hitPhi >= (int)omtfConfig->nPhiBins())
        continue;

      auto& hitTree = layerTree.add("Hit", "");

      hitTree.add("<xmlattr>.iInput", iHit);
      hitTree.add("<xmlattr>.iEta", input->getHitEta(iLayer, iHit));
      hitTree.add("<xmlattr>.iPhi", hitPhi);

      //in the firmware the hit quality is taken from link data only for the DT stubs
      //for the CSC and RPC 1 means the hit is valid, 0 - not.
      //in the input it is still worth to have the actual quality of the CSC and RPC
      //Because it might be used in the neural network
      if (iLayer >= 6)
        hitTree.add("<xmlattr>.iQual", 1);
      else
        hitTree.add("<xmlattr>.iQual", input->getHitQual(iLayer, iHit));
    }

    if (!layerTree.empty()) {
      layerTree.add("<xmlattr>.iLayer", iLayer);
      procTree.add_child("Layer", layerTree);
    }
  }
  */

  for (auto& algoCand : algoCandidates) {
    ///Dump only regions, where a candidate was found
    if (algoCand->isValid()) {
      auto& algoMuonTree = procTree.add("AlgoMuon", "");
      algoMuonTree.add("<xmlattr>.charge", algoCand->getChargeConstr());
      algoMuonTree.add("<xmlattr>.disc", algoCand->getDisc());
      algoMuonTree.add("<xmlattr>.pdfSumConstr", algoCand->getGpResultConstr().getPdfSum());
      algoMuonTree.add("<xmlattr>.pdfSumUnconstr", algoCand->getGpResultUnconstr().getPdfSumUnconstr());
      algoMuonTree.add("<xmlattr>.etaCode", algoCand->getEtaHw());
      algoMuonTree.add("<xmlattr>.iRefHit", algoCand->getRefHitNumber());
      algoMuonTree.add("<xmlattr>.iRefLayer", algoCand->getRefLayer());

      //algoMuonTree.add("<xmlattr>.layers", std::bitset<18>(algoCand->getFiredLayerBits()));

      algoMuonTree.add("<xmlattr>.layersConstr", std::bitset<18>(algoCand->getGpResultConstr().getFiredLayerBits()));
      algoMuonTree.add("<xmlattr>.layersUnconstr",
                       std::bitset<18>(algoCand->getGpResultUnconstr().getFiredLayerBits()));

      algoMuonTree.add("<xmlattr>.nHits", algoCand->getQ());

      algoMuonTree.add("<xmlattr>.firedCntConstr", algoCand->getGpResultConstr().getFiredLayerCnt());
      algoMuonTree.add("<xmlattr>.firedCntUnconstr", algoCand->getGpResultUnconstr().getFiredLayerCnt());

      algoMuonTree.add("<xmlattr>.patNumConstr", algoCand->getHwPatternNumConstr());
      algoMuonTree.add("<xmlattr>.patNumUnconstr", algoCand->getHwPatternNumUnconstr());

      algoMuonTree.add("<xmlattr>.phiCode", algoCand->getPhi());

      algoMuonTree.add("<xmlattr>.phiConstr", algoCand->getGpResultConstr().getPhi());
      algoMuonTree.add("<xmlattr>.phiUnConstr", algoCand->getGpResultUnconstr().getPhi());

      algoMuonTree.add("<xmlattr>.phiRHit", algoCand->getPhiRHit());

      //in the firmware, the algoMuon has no pt nor upt yet,
      //only the pattern number, which is converted to the hwpt in the ghostbuster
      algoMuonTree.add("<xmlattr>.ptCodeConstr", algoCand->getPtConstr());
      algoMuonTree.add("<xmlattr>.ptCodeUnconstr", algoCand->getPtUnconstr());

      auto& gpResultTree = algoMuonTree.add("gpResultConstr", "");
      auto& gpResultConstr = algoCand->getGpResultConstr();

      gpResultTree.add("<xmlattr>.patNum", algoCand->getHwPatternNumConstr());
      gpResultTree.add("<xmlattr>.pdfSum", gpResultConstr.getPdfSum());

      for (unsigned int iLogicLayer = 0; iLogicLayer < gpResultConstr.getStubResults().size(); ++iLogicLayer) {
        auto& layerTree = gpResultTree.add("layer", "");
        layerTree.add("<xmlattr>.num", iLogicLayer);
        auto pdfBin = gpResultConstr.getStubResults()[iLogicLayer].getPdfBin();
        if (pdfBin == 5400)
          pdfBin = 0;
        layerTree.add("<xmlattr>.pdfBin", pdfBin);
        layerTree.add("<xmlattr>.pdfVal", gpResultConstr.getStubResults()[iLogicLayer].getPdfVal());
        layerTree.add("<xmlattr>.fired", gpResultConstr.isLayerFired(iLogicLayer));
      }

      if (algoCand->getGpResultUnconstr().isValid()) {
        auto& gpResultTree = algoMuonTree.add("gpResultUnconstr", "");
        auto& gpResult = algoCand->getGpResultUnconstr();

        gpResultTree.add("<xmlattr>.patNum", algoCand->getHwPatternNumUnconstr());
        gpResultTree.add("<xmlattr>.pdfSum", gpResult.getPdfSumUnconstr());

        for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
          auto& layerTree = gpResultTree.add("layer", "");
          layerTree.add("<xmlattr>.num", iLogicLayer);
          auto pdfBin = gpResult.getStubResults()[iLogicLayer].getPdfBin();
          if (pdfBin == 5400)
            pdfBin = 0;
          layerTree.add("<xmlattr>.pdfBin", pdfBin);
          layerTree.add("<xmlattr>.pdfVal", gpResult.getStubResults()[iLogicLayer].getPdfVal());
          layerTree.add("<xmlattr>.fired", gpResult.isLayerFired(iLogicLayer));
        }
      }
    }
  }

  // WP1A Phase C (2026-09-14): GhostBuster output -- algoCandidates above is
  // GhostBuster's real INPUT (pattern-match results, one per fired refHit,
  // pre-ghost-killing); gbCandidates is its real OUTPUT (survivors after
  // GhostBuster::select()/killGhosts(), same AlgoMuon type, same real
  // fields) and had no XML representation at all before this -- the "OMTF
  // GhostBuster in/out" row action-plan section 5.4 flagged as "new".
  // Same exact field set as AlgoMuon above (kept duplicated, not factored
  // into a shared helper, matching this codebase's own established
  // barrel/endcap-style duplication pattern) so a real ghost_buster HLS
  // module trace can be compared directly against this real ground truth.
  for (auto& gbCand : gbCandidates) {
    if (gbCand->isValid()) {
      auto& gbMuonTree = procTree.add("GbMuon", "");
      gbMuonTree.add("<xmlattr>.charge", gbCand->getChargeConstr());
      gbMuonTree.add("<xmlattr>.disc", gbCand->getDisc());
      gbMuonTree.add("<xmlattr>.pdfSumConstr", gbCand->getGpResultConstr().getPdfSum());
      gbMuonTree.add("<xmlattr>.pdfSumUnconstr", gbCand->getGpResultUnconstr().getPdfSumUnconstr());
      gbMuonTree.add("<xmlattr>.etaCode", gbCand->getEtaHw());
      gbMuonTree.add("<xmlattr>.iRefHit", gbCand->getRefHitNumber());
      gbMuonTree.add("<xmlattr>.iRefLayer", gbCand->getRefLayer());

      gbMuonTree.add("<xmlattr>.layersConstr", std::bitset<18>(gbCand->getGpResultConstr().getFiredLayerBits()));
      gbMuonTree.add("<xmlattr>.layersUnconstr",
                     std::bitset<18>(gbCand->getGpResultUnconstr().getFiredLayerBits()));

      gbMuonTree.add("<xmlattr>.nHits", gbCand->getQ());

      gbMuonTree.add("<xmlattr>.firedCntConstr", gbCand->getGpResultConstr().getFiredLayerCnt());
      gbMuonTree.add("<xmlattr>.firedCntUnconstr", gbCand->getGpResultUnconstr().getFiredLayerCnt());

      gbMuonTree.add("<xmlattr>.patNumConstr", gbCand->getHwPatternNumConstr());
      gbMuonTree.add("<xmlattr>.patNumUnconstr", gbCand->getHwPatternNumUnconstr());

      gbMuonTree.add("<xmlattr>.phiCode", gbCand->getPhi());

      gbMuonTree.add("<xmlattr>.phiConstr", gbCand->getGpResultConstr().getPhi());
      gbMuonTree.add("<xmlattr>.phiUnConstr", gbCand->getGpResultUnconstr().getPhi());

      gbMuonTree.add("<xmlattr>.phiRHit", gbCand->getPhiRHit());

      gbMuonTree.add("<xmlattr>.ptCodeConstr", gbCand->getPtConstr());
      gbMuonTree.add("<xmlattr>.ptCodeUnconstr", gbCand->getPtUnconstr());
    }
  }

  for (auto& finalMuon : finalMuons) {
    auto& candMuonTree = procTree.add("CandMuon", "");
    candMuonTree.add("<xmlattr>.hwEta", finalMuon->getEtaGmt());
    candMuonTree.add("<xmlattr>.hwPhi", finalMuon->getPhiGmt());
    candMuonTree.add("<xmlattr>.hwPt", finalMuon->getPtGmt());
    candMuonTree.add("<xmlattr>.hwUPt", finalMuon->getPtUnconstrGmt());
    candMuonTree.add("<xmlattr>.hwQual", finalMuon->getQuality());
    candMuonTree.add("<xmlattr>.hwSign", finalMuon->getSign());
    candMuonTree.add("<xmlattr>.hwSignValid", 1);
    candMuonTree.add("<xmlattr>.hwTrackAddress", std::bitset<29>(finalMuon->getAlgoMuon()->getFiredLayerBits()));
    candMuonTree.add("<xmlattr>.link", (mtfType == l1t::omtf_neg ? 60 + iProcessor : 42 + iProcessor));
    candMuonTree.add("<xmlattr>.processor", iProcessor);

    std::ostringstream stringStr;
    if (mtfType == l1t::omtf_neg)
      stringStr << "OMTF_NEG";
    else if (mtfType == l1t::omtf_pos)
      stringStr << "OMTF_POS";

    candMuonTree.add("<xmlattr>.trackFinderType", stringStr.str());
  }

  // WP1A Phase C (2026-09-14): SAMuon-equivalent export -- replicates
  // OmtfProcessorPhase2::getSAMuons()'s real selection/fallback logic
  // exactly (same iProcessor/mtfType/finalMuons granularity that function
  // is itself called at in OmtfEmulation::run()), producing the same two
  // aligned constrained/unconstrained collections real GMT-facing SAMuon
  // objects would form -- something CandMuon's single list doesn't
  // capture (CandMuon predates the Phase-2 constrained/unconstrained
  // SAMuon split and has no equivalent of getSAMuons()'s "backfill the
  // unconstrained collection with the constrained candidate when no real
  // unconstrained measurement exists, so both collections stay the same
  // size" rule).
  //
  // IMPORTANT, confirmed by reading getSAMuons()/SAMuon.h directly, not
  // assumed: hwD0/hwZ0/hwBeta/word() are NOT real per-candidate
  // measurements in CMSSW_20_0_0 today.
  //  - z0 is unconditionally 0 (no computation exists anywhere).
  //  - d0 is a fixed placeholder CMSSW's own comment states outright:
  //    "Set d0 to the default values of zero (0) and a large number (50)
  //    ... until the d0 measurement is implemented" -- 0 for constrained,
  //    50/Phase2L1GMT::LSBSAd0 (=13 in hw units, LSBSAd0=3.84) for
  //    unconstrained, replicated exactly below (not re-derived).
  //  - word() (the real packed 64-bit GMT word) is dead code in
  //    OmtfProcessorPhase2.cc -- the whole packing block and its
  //    setWord() call are commented out.
  //  - hwBeta has no setter called anywhere in getSAMuons(); the SAMuon
  //    constructor used here doesn't even take it.
  // So this export deliberately does NOT include hwD0/hwZ0/hwBeta/word --
  // exporting a fixed constant or an always-zero field as a "golden
  // reference" would misrepresent CMSSW's own not-yet-implemented state as
  // if it were real ground truth. pt/eta/phi/charge/quality ARE real
  // (same fields CandMuon already carries); what's genuinely new here is
  // WHICH candidates appear in each collection and with what pt.
  for (auto& finalMuon : finalMuons) {
    int ptConstr = finalMuon->getPtGmt();
    if (ptConstr > 0) {
      auto& saMuonTree = procTree.add("SAMuonConstr", "");
      saMuonTree.add("<xmlattr>.hwPt", ptConstr);
      saMuonTree.add("<xmlattr>.hwEta", finalMuon->getEtaGmt());
      saMuonTree.add("<xmlattr>.hwPhi", finalMuon->getPhiGmt());
      saMuonTree.add("<xmlattr>.hwSign", finalMuon->getSign());
      saMuonTree.add("<xmlattr>.hwQual", finalMuon->getQuality());
    }

    int ptUnconstr = finalMuon->getPtUnconstrGmt();
    if (ptUnconstr == 0) {
      // real getSAMuons() fallback: no unconstrained measurement -> use
      // the constrained candidate instead, so both real SAMuon
      // collections stay the same size (their own comment's own words).
      ptUnconstr = ptConstr;
    }
    if (ptUnconstr > 0) {
      auto& saMuonTree = procTree.add("SAMuonUnconstr", "");
      saMuonTree.add("<xmlattr>.hwPt", ptUnconstr);
      saMuonTree.add("<xmlattr>.hwEta", finalMuon->getEtaGmt());
      saMuonTree.add("<xmlattr>.hwPhi", finalMuon->getPhiGmt());
      saMuonTree.add("<xmlattr>.hwSign", finalMuon->getSign());
      saMuonTree.add("<xmlattr>.hwQual", finalMuon->getQuality());
    }
  }

  if (!procTree.empty())
    eventTree->add_child("Processor", procTree);
}

void XMLEventWriter::observeEventBegin(const edm::Event& iEvent) {
  eventNum++;
  eventId = iEvent.id().event();

  eventTree = &(tree.add("OMTF.Event", ""));
  eventTree->add("<xmlattr>.iEvent", eventId);

  eventTree = &(eventTree->add("bx", ""));
  eventTree->add("<xmlattr>.iBx", 2 * eventId);
}

void XMLEventWriter::observeEventEnd(const edm::Event& iEvent, FinalMuons& finalMuons) {
  if (eventsPerFile > 0 && eventNum > 0 && (eventNum % (unsigned int)eventsPerFile) == 0) {
    flushCurrentTreeToFile();
  }
}

void XMLEventWriter::endJob() {
  edm::LogInfo("l1tOmtfEventPrint") << "XMLEventWriter::endJob() - writing the data to the xml - starting";

  std::string outputName = fName;
  if (eventsPerFile > 0) {
    // Write the final (possibly partial) batch with a part index as well
    std::string baseName = fName;
    std::string ext;
    auto dotPos = fName.rfind('.');
    if (dotPos != std::string::npos) {
      baseName = fName.substr(0, dotPos);
      ext = fName.substr(dotPos);
    }
    std::ostringstream ss;
    ss << baseName << "_part" << std::setfill('0') << std::setw(3) << fileIndex << ext;
    outputName = ss.str();
  }

  boost::property_tree::write_xml(
      outputName, tree, std::locale(), boost::property_tree::xml_parser::xml_writer_make_settings<std::string>(' ', 2));
  edm::LogInfo("l1tOmtfEventPrint") << "XMLEventWriter::endJob() - writing the data to the xml - done";
}
