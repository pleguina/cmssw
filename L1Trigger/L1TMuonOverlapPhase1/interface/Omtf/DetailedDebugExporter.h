/*
 * DetailedDebugExporter.h
 *
 * Exports detailed debug information for a specific event,
 * capturing all StubResult values from process1Layer1RefLayer calls
 * and finalise function results for deep algorithm analysis.
 */

#ifndef L1T_OmtfP1_DETAILEDDEBUGEXPORTER_H_
#define L1T_OmtfP1_DETAILEDDEBUGEXPORTER_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/StubResult.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternResult.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/FinalMuon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/multi_array.hpp>
#include <map>
#include <string>
#include <sstream>

class DetailedDebugExporter : public IOMTFEmulationObserver {
public:
  DetailedDebugExporter(const OMTFConfiguration* omtfConfig,
                        int targetEventNumber,
                        const std::string& outputDir = "./")
      : omtfConfig_(omtfConfig),
        targetEventNumber_(targetEventNumber),
        outputDir_(outputDir),
        currentEventNumber_(-1),
        isTargetEvent_(false) {}

  ~DetailedDebugExporter() override = default;

  // Required pure virtual methods from base class
  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>& input,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const FinalMuons& finalMuons) override {
    // Not used for detailed debug
  }

  void endJob() override {
    // Nothing to do at end of job
  }

  void observeEventBegin(const edm::Event& iEvent) override {
    currentEventNumber_ = iEvent.id().event();
    isTargetEvent_ = (currentEventNumber_ == targetEventNumber_);

    if (isTargetEvent_) {
      edm::LogInfo("DetailedDebugExporter") << "Target event " << targetEventNumber_
                                            << " detected - capturing detailed debug info";
      debugTree_ = boost::property_tree::ptree();
      debugTree_.add("<xmlattr>.eventNumber", currentEventNumber_);
      debugTree_.add("<xmlattr>.run", iEvent.id().run());
      debugTree_.add("<xmlattr>.lumi", iEvent.id().luminosityBlock());
    }
  }

  void observeEventEnd(const edm::Event& iEvent,
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override {
    if (isTargetEvent_) {
      // Write the detailed debug XML
      writeDebugXML();
      edm::LogInfo("DetailedDebugExporter") << "Detailed debug XML written for event " << currentEventNumber_;
    }
  }

  // Capture process1Layer1RefLayer results
  void observeStubResult(unsigned int iProcessor,
                         unsigned int iRefHit,
                         unsigned int iRefLayer,
                         unsigned int iLayer,
                         unsigned int patternNumber,
                         const Key& patternKey,
                         const StubResult& stubResult,
                         const MuonStubPtrs1D& layerStubs,
                         const std::vector<int>& extrapolatedPhi,
                         const MuonStubPtr& refStub) {
    if (!isTargetEvent_)
      return;

    // Build hierarchical path: processor -> refHit -> pattern -> layer
    std::ostringstream pathKey;
    pathKey << "proc" << iProcessor << "_refHit" << iRefHit << "_pattern" << patternNumber;

    boost::property_tree::ptree stubResultNode;

    // Add pattern info (only once per pattern)
    if (layerResults_[pathKey.str()].empty()) {
      boost::property_tree::ptree& processorNode = getOrCreateNode("processor", iProcessor);
      boost::property_tree::ptree& refHitNode = getOrCreateNode(processorNode, "refHit", iRefHit);
      boost::property_tree::ptree& patternNode = getOrCreateNode(refHitNode, "pattern", patternNumber);

      patternNode.add("<xmlattr>.patternNumber", patternNumber);
      patternNode.add("<xmlattr>.pt", patternKey.thePt);
      patternNode.add("<xmlattr>.charge", patternKey.theCharge);
      patternNode.add("<xmlattr>.etaCode", patternKey.theEtaCode);
      patternNode.add("<xmlattr>.iRefLayer", iRefLayer);

      // Add reference stub info
      if (refStub) {
        boost::property_tree::ptree refStubNode;
        refStubNode.add("<xmlattr>.phiHw", refStub->phiHw);
        refStubNode.add("<xmlattr>.etaHw", refStub->etaHw);
        refStubNode.add("<xmlattr>.phiBHw", refStub->phiBHw);
        refStubNode.add("<xmlattr>.quality", refStub->qualityHw);
        refStubNode.add("<xmlattr>.logicLayer", refStub->logicLayer);
        refStubNode.add("<xmlattr>.detId", refStub->detId);
        patternNode.add_child("refStub", refStubNode);
      }
    }

    // Add layer result
    stubResultNode.add("<xmlattr>.layer", iLayer);
    stubResultNode.add("<xmlattr>.valid", stubResult.getValid());
    stubResultNode.add("<xmlattr>.pdfVal", stubResult.getPdfVal());
    stubResultNode.add("<xmlattr>.pdfBin", stubResult.getPdfBin());
    stubResultNode.add("<xmlattr>.deltaPhi", stubResult.getDeltaPhi());

    // Add stub info if present
    if (stubResult.getMuonStub()) {
      const auto& stub = stubResult.getMuonStub();
      stubResultNode.add("<xmlattr>.stubPhiHw", stub->phiHw);
      stubResultNode.add("<xmlattr>.stubEtaHw", stub->etaHw);
      stubResultNode.add("<xmlattr>.stubPhiBHw", stub->phiBHw);
      stubResultNode.add("<xmlattr>.stubQuality", stub->qualityHw);
      stubResultNode.add("<xmlattr>.stubLogicLayer", stub->logicLayer);
    }

    // Add extrapolated phi if available
    if (!extrapolatedPhi.empty() && stubResult.getMuonStub()) {
      // Find the stub index in layerStubs
      for (size_t iStub = 0; iStub < layerStubs.size(); ++iStub) {
        if (layerStubs[iStub] && layerStubs[iStub] == stubResult.getMuonStub()) {
          if (iStub < extrapolatedPhi.size()) {
            stubResultNode.add("<xmlattr>.extrapolatedPhi", extrapolatedPhi[iStub]);
          }
          break;
        }
      }
    }

    layerResults_[pathKey.str()].add_child("layerResult", stubResultNode);
  }

  // Capture finalise results - takes GoldenPatternResult array instead of AlgoMuons
  void observeFinaliseResults(unsigned int iProcessor,
                              unsigned int patternNumber,
                              const Key& patternKey,
                              const boost::detail::multi_array::sub_array<GoldenPatternResult, 1>& gpResults) {
    if (!isTargetEvent_)
      return;

    boost::property_tree::ptree& processorNode = getOrCreateNode("processor", iProcessor);
    boost::property_tree::ptree finaliseNode;

    finaliseNode.add("<xmlattr>.patternNumber", patternNumber);
    finaliseNode.add("<xmlattr>.pt", patternKey.thePt);
    finaliseNode.add("<xmlattr>.charge", patternKey.theCharge);
    finaliseNode.add("<xmlattr>.etaCode", patternKey.theEtaCode);

    // Add GoldenPatternResult info for each refHit
    for (size_t iRefHit = 0; iRefHit < gpResults.size(); ++iRefHit) {
      const auto& gpResult = gpResults[iRefHit];
      if (gpResult.isValid()) {
        boost::property_tree::ptree resultNode;
        resultNode.add("<xmlattr>.refHit", iRefHit);
        resultNode.add("<xmlattr>.refLayer", gpResult.getRefLayer());
        resultNode.add("<xmlattr>.phi", gpResult.getPhi());
        resultNode.add("<xmlattr>.eta", gpResult.getEta());
        resultNode.add("<xmlattr>.refHitPhi", gpResult.getRefHitPhi());
        resultNode.add("<xmlattr>.pdfSum", gpResult.getPdfSum());
        resultNode.add("<xmlattr>.pdfSumUnconstr", gpResult.getPdfSumUnconstr());
        resultNode.add("<xmlattr>.firedLayerCnt", gpResult.getFiredLayerCnt());
        resultNode.add("<xmlattr>.firedLayerBits", gpResult.getFiredLayerBits());
        resultNode.add("<xmlattr>.gpProbability1", gpResult.getGpProbability1());
        resultNode.add("<xmlattr>.gpProbability2", gpResult.getGpProbability2());

        finaliseNode.add_child("gpResult", resultNode);
      }
    }

    processorNode.add_child("finaliseResult", finaliseNode);
  }

private:
  boost::property_tree::ptree& getOrCreateNode(const std::string& nodeType, unsigned int index) {
    std::ostringstream key;
    key << nodeType << index;

    // Check if node already exists
    for (auto& child : debugTree_) {
      if (child.first == nodeType) {
        auto indexAttr = child.second.get_optional<unsigned int>("<xmlattr>.index");
        if (indexAttr && *indexAttr == index) {
          return child.second;
        }
      }
    }

    // Create new node
    boost::property_tree::ptree newNode;
    newNode.add("<xmlattr>.index", index);
    return debugTree_.add_child(nodeType, newNode);
  }

  boost::property_tree::ptree& getOrCreateNode(boost::property_tree::ptree& parent,
                                                 const std::string& nodeType,
                                                 unsigned int index) {
    // Check if node already exists
    for (auto& child : parent) {
      if (child.first == nodeType) {
        auto indexAttr = child.second.get_optional<unsigned int>("<xmlattr>.index");
        if (indexAttr && *indexAttr == index) {
          return child.second;
        }
      }
    }

    // Create new node
    boost::property_tree::ptree newNode;
    newNode.add("<xmlattr>.index", index);
    return parent.add_child(nodeType, newNode);
  }

  void writeDebugXML() {
    // Merge layer results into debug tree
    for (const auto& entry : layerResults_) {
      // Parse the key: procX_refHitY_patternZ
      std::string key = entry.first;
      size_t procPos = key.find("proc");
      size_t refHitPos = key.find("_refHit");
      size_t patternPos = key.find("_pattern");

      if (procPos == std::string::npos || refHitPos == std::string::npos || patternPos == std::string::npos)
        continue;

      unsigned int iProc = std::stoi(key.substr(procPos + 4, refHitPos - (procPos + 4)));
      unsigned int iRefHit = std::stoi(key.substr(refHitPos + 7, patternPos - (refHitPos + 7)));
      unsigned int iPattern = std::stoi(key.substr(patternPos + 8));

      boost::property_tree::ptree& procNode = getOrCreateNode("processor", iProc);
      boost::property_tree::ptree& refHitNode = getOrCreateNode(procNode, "refHit", iRefHit);
      boost::property_tree::ptree& patternNode = getOrCreateNode(refHitNode, "pattern", iPattern);

      // Add all layer results to this pattern
      for (const auto& layerResult : entry.second) {
        patternNode.add_child(layerResult.first, layerResult.second);
      }
    }

    // Write to file
    std::ostringstream filename;
    filename << outputDir_;
    if (outputDir_.back() != '/')
      filename << "/";
    filename << "DetailedDebug_Event" << currentEventNumber_ << ".xml";

    boost::property_tree::ptree eventTree;
    eventTree.add_child("Event", debugTree_);

    auto settings = boost::property_tree::xml_writer_make_settings<std::string>('\t', 1);
    boost::property_tree::write_xml(filename.str(), eventTree, std::locale(), settings);

    edm::LogInfo("DetailedDebugExporter") << "Detailed debug XML written to: " << filename.str();
  }

  const OMTFConfiguration* omtfConfig_;
  int targetEventNumber_;
  std::string outputDir_;
  int currentEventNumber_;
  bool isTargetEvent_;

  boost::property_tree::ptree debugTree_;
  std::map<std::string, boost::property_tree::ptree> layerResults_;  // Temporary storage for layer results
};

#endif /* L1T_OmtfP1_DETAILEDDEBUGEXPORTER_H_ */
