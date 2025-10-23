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
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternBase.h"
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

    // Match based on event ID (not sequential counter)
    isTargetEvent_ = (currentEventNumber_ == targetEventNumber_);

    if (isTargetEvent_) {
      edm::LogInfo("DetailedDebugExporter") << "Target event ID=" << currentEventNumber_
                                            << " detected - capturing detailed debug info";
      debugTree_ = boost::property_tree::ptree();
      // Don't add attributes here - will be added in writeDebugXML with proper structure
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

  // Capture restricted stubs for each refHit
  void observeRestrictedStubs(unsigned int iProcessor,
                             unsigned int iRefHit,
                             unsigned int iRefLayer,
                             const MuonStubPtrs2D& restrictedStubs,
                             const std::vector<std::pair<unsigned int, std::vector<int>>>& extrapolatedPhis,
                             const MuonStubPtr& refStub) {
    if (!isTargetEvent_)
      return;

    // Build hierarchical path for this refHit
    std::ostringstream refHitKey;
    refHitKey << "proc" << iProcessor << "_refHit" << iRefHit;

    // Store restricted stubs for this refHit (will be added to XML later)
    restrictedStubsPerRefHit_[refHitKey.str()] = std::make_tuple(iRefLayer, restrictedStubs, extrapolatedPhis, refStub);
    
    std::cout << "DetailedDebugExporter: Stored restricted stubs for " << refHitKey.str() 
              << " iRefLayer=" << iRefLayer 
              << " nLayers=" << restrictedStubs.size() << std::endl;
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
                         const MuonStubPtr& refStub,
                         const GoldenPatternBase* goldenPattern) {
    if (!isTargetEvent_)
      return;

    // Build hierarchical path: processor -> refHit -> pattern -> layer
    std::ostringstream pathKey;
    pathKey << "proc" << iProcessor << "_refHit" << iRefHit << "_pattern" << patternNumber;

    boost::property_tree::ptree stubResultNode;

    // Add pattern info (only once per pattern)
    if (layerResults_[pathKey.str()].empty()) {
      boost::property_tree::ptree& processorNode = getOrCreateNode("processor", iProcessor);
      boost::property_tree::ptree& refHitNode = getOrCreateNode(processorNode, "referenceHit", iRefHit);
      boost::property_tree::ptree& patternNode = getOrCreateNode(refHitNode, "pattern", patternNumber);

      patternNode.add("<xmlattr>.patternNumber", patternNumber);
      patternNode.add("<xmlattr>.pt", patternKey.thePt);
      patternNode.add("<xmlattr>.charge", patternKey.theCharge);
      patternNode.add("<xmlattr>.etaCode", patternKey.theEtaCode);
      patternNode.add("<xmlattr>.iRefLayer", iRefLayer);

      // Add reference stub info
      if (refStub) {
        boost::property_tree::ptree refStubNode;
        refStubNode.add("<xmlattr>.refPhi", refStub->phiHw);
        refStubNode.add("<xmlattr>.refEta", refStub->etaHw);
        refStubNode.add("<xmlattr>.refPhiB", refStub->phiBHw);
        refStubNode.add("<xmlattr>.refQuality", refStub->qualityHw);
        refStubNode.add("<xmlattr>.refLogicLayer", refStub->logicLayer);
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
    
    // Add phiDistMin (raw distance before pdfMiddle offset)
    // pdfMiddle = 1 << (nPdfAddrBits - 1) = 64 for nPdfAddrBits=7
    int pdfMiddle = 1 << (omtfConfig_->nPdfAddrBits() - 1);
    int phiDistMin = stubResult.getPdfBin() - pdfMiddle;
    stubResultNode.add("<xmlattr>.phiDistMin", phiDistMin);

    // Add meanDistPhi and shift from the golden pattern
    if (goldenPattern && refStub) {
      int meanDistPhi = goldenPattern->meanDistPhiValue(iLayer, iRefLayer, refStub->phiBHw);
      int shift = goldenPattern->getDistPhiBitShift(iLayer, iRefLayer);
      stubResultNode.add("<xmlattr>.meanDistPhi", meanDistPhi);
      stubResultNode.add("<xmlattr>.shift", shift);
    }

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
    // First, add restricted stubs to refHit nodes
    for (const auto& entry : restrictedStubsPerRefHit_) {
      // Parse the key: procX_refHitY
      std::string key = entry.first;
      size_t procPos = key.find("proc");
      size_t refHitPos = key.find("_refHit");

      if (procPos == std::string::npos || refHitPos == std::string::npos)
        continue;

      unsigned int iProc = std::stoi(key.substr(procPos + 4, refHitPos - (procPos + 4)));
      unsigned int iRefHit = std::stoi(key.substr(refHitPos + 7));

      const auto& data = entry.second;
      unsigned int iRefLayer = std::get<0>(data);
      const auto& restrictedStubs = std::get<1>(data);
      const auto& extrapolatedPhis = std::get<2>(data);
      const auto& refStub = std::get<3>(data);

      boost::property_tree::ptree& procNode = getOrCreateNode("processor", iProc);
      boost::property_tree::ptree& refHitNode = getOrCreateNode(procNode, "referenceHit", iRefHit);

      // Add reference stub attributes directly to referenceHit tag
      refHitNode.add("<xmlattr>.iRefLayer", iRefLayer);
      if (refStub) {
        refHitNode.add("<xmlattr>.refPhi", refStub->phiHw);
        refHitNode.add("<xmlattr>.refEta", refStub->etaHw);
        refHitNode.add("<xmlattr>.refPhiB", refStub->phiBHw);
        refHitNode.add("<xmlattr>.refQuality", refStub->qualityHw);
        refHitNode.add("<xmlattr>.refLogicLayer", refStub->logicLayer);
      }

      // Add restricted stubs section
      boost::property_tree::ptree restrictedStubsNode;

      // Add stubs from all layers
      for (unsigned int iLayer = 0; iLayer < restrictedStubs.size(); ++iLayer) {
        const auto& layerStubs = restrictedStubs[iLayer];
        
        // Find extrapolated phi for this layer
        std::vector<int> layerExtrapolatedPhi;
        for (const auto& pair : extrapolatedPhis) {
          if (pair.first == iLayer) {
            layerExtrapolatedPhi = pair.second;
            break;
          }
        }

        for (size_t iStub = 0; iStub < layerStubs.size(); ++iStub) {
          const auto& stub = layerStubs[iStub];
          if (stub) {
            boost::property_tree::ptree stubNode;
            stubNode.add("<xmlattr>.iLayer", iLayer);
            stubNode.add("<xmlattr>.inputNumber", iStub);
            stubNode.add("<xmlattr>.phi", stub->phiHw);
            stubNode.add("<xmlattr>.phiB", stub->phiBHw);
            stubNode.add("<xmlattr>.eta", stub->etaHw);
            stubNode.add("<xmlattr>.quality", stub->qualityHw);
            stubNode.add("<xmlattr>.logicLayer", stub->logicLayer);
            stubNode.add("<xmlattr>.detId", stub->detId);
            
            if (iStub < layerExtrapolatedPhi.size()) {
              stubNode.add("<xmlattr>.extrapolatedPhi", layerExtrapolatedPhi[iStub]);
            }
            
            restrictedStubsNode.add_child("stub", stubNode);
          }
        }
      }

      refHitNode.add_child("restrictedStubs", restrictedStubsNode);
    }

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
      boost::property_tree::ptree& refHitNode = getOrCreateNode(procNode, "referenceHit", iRefHit);
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
    filename << "DetailedSingleEv_" << currentEventNumber_ << ".xml";

    // Create XML structure matching TestEvents.xml:
    // <OMTF version="0x0000">
    //   <Event iEvent="X">
    //     <bx iBx="0">
    //       <Processor board="OMTFpX" iProcessor="X" position="+/-1">
    //         ... processor data ...
    
    // Wrap each processor in its own Processor node with bx
    boost::property_tree::ptree eventTree;
    eventTree.add("<xmlattr>.iEvent", currentEventNumber_);
    
    // For detailed debug, assume bx=0 (central BX)
    boost::property_tree::ptree bxTree;
    bxTree.add("<xmlattr>.iBx", 0);
    
    // Move all processor nodes into Processor tags with proper attributes
    for (auto& procNode : debugTree_) {
      if (procNode.first == "processor") {
        unsigned int iProc = procNode.second.get<unsigned int>("<xmlattr>.index");
        
        boost::property_tree::ptree processorTree;
        std::ostringstream boardName;
        boardName << "OMTFp" << iProc;
        processorTree.add("<xmlattr>.board", boardName.str());
        processorTree.add("<xmlattr>.iProcessor", iProc);
        processorTree.add("<xmlattr>.position", "+1");  // Default to positive, adjust if needed
        
        // Create referenceHits wrapper
        boost::property_tree::ptree referenceHitsTree;
        
        // Move all children from the processor node
        for (auto& child : procNode.second) {
          if (child.first != "<xmlattr>") {
            // Wrap referenceHit elements in referenceHits container
            if (child.first == "referenceHit") {
              referenceHitsTree.add_child("referenceHit", child.second);
            } else {
              processorTree.add_child(child.first, child.second);
            }
          }
        }
        
        // Add the referenceHits wrapper if it has children
        if (!referenceHitsTree.empty()) {
          processorTree.add_child("referenceHits", referenceHitsTree);
        }
        
        bxTree.add_child("Processor", processorTree);
      }
    }
    
    eventTree.add_child("bx", bxTree);
    
    // Create OMTF element with Event child
    boost::property_tree::ptree omtfElement;
    omtfElement.add("<xmlattr>.version", "0x0000");
    omtfElement.add_child("Event", eventTree);
    
    // Create root tree with OMTF as child
    boost::property_tree::ptree rootTree;
    rootTree.add_child("OMTF", omtfElement);
    
    auto settings = boost::property_tree::xml_writer_make_settings<std::string>('\t', 1);
    boost::property_tree::write_xml(filename.str(), rootTree, std::locale(), settings);

    edm::LogInfo("DetailedDebugExporter") << "Detailed debug XML written to: " << filename.str();
  }

  const OMTFConfiguration* omtfConfig_;
  int targetEventNumber_;  // Target event ID to capture
  std::string outputDir_;
  int currentEventNumber_;  // Current event ID from file
  bool isTargetEvent_;

  boost::property_tree::ptree debugTree_;
  std::map<std::string, boost::property_tree::ptree> layerResults_;  // Temporary storage for layer results
  
  // Storage for restricted stubs per refHit: key="procX_refHitY", value=tuple(iRefLayer, stubs, extrapolatedPhis, refStub)
  std::map<std::string, std::tuple<unsigned int, MuonStubPtrs2D, std::vector<std::pair<unsigned int, std::vector<int>>>, MuonStubPtr>> restrictedStubsPerRefHit_;
};

#endif /* L1T_OmtfP1_DETAILEDDEBUGEXPORTER_H_ */
