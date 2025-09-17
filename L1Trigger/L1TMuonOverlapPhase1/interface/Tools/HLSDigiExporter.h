#ifndef L1T_OmtfP1_HLSDIGIEXPORTER_H_
#define L1T_OmtfP1_HLSDIGIEXPORTER_H_

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/IOMTFEmulationObserver.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/MuonStub.h"
#include "FWCore/Framework/interface/Event.h"

#include <fstream>
#include <string>
#include <memory>
#include <vector>

/**
 * HLS Digi Exporter
 * 
 * This class captures ALL input digis (DT, CSC, RPC) and exports them to CSV files
 * for HLS test vector generation. It works alongside the existing XML system.
 */

class HLSDigiExporter : public IOMTFEmulationObserver {
public:
  HLSDigiExporter(const std::string& outputDir = "hls_input_digis");
  ~HLSDigiExporter() override;

  // Observer interface
  void observeEventBegin(const edm::Event& iEvent) override;
  void observeEventEnd(const edm::Event& iEvent, 
                       std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) override;
  
  void observeProcesorBegin(unsigned int iProcessor, l1t::tftype mtfType) override;
  
  void addProcesorData(std::string key, boost::property_tree::ptree& procDataTree) override;
  
  void observeProcesorEmulation(unsigned int iProcessor,
                                l1t::tftype mtfType,
                                const std::shared_ptr<OMTFinput>& input,
                                const AlgoMuons& algoCandidates,
                                const AlgoMuons& gbCandidates,
                                const FinalMuons& finalMuons) override {}

  // New method to observe reference hit processing
  void observeRefHitProcessing(unsigned int iProcessor,
                               unsigned int iRefHit,
                               const RefHitDef& refHitDef,
                               const std::vector<std::pair<unsigned int, MuonStubPtrs1D>>& allLayerStubs);

  void endJob() override;

private:
  std::string outputDir_;
  unsigned int currentEvent_;
  unsigned int currentRun_;
  unsigned int currentProcessor_;
  l1t::tftype currentMtfType_;
  
  // CSV files for different digi types
  std::ofstream dtPhiDigiFile_;
  std::ofstream dtThetaDigiFile_;
  std::ofstream cscDigiFile_;
  std::ofstream rpcDigiFile_;
  std::ofstream goldenResultsFile_;  // For converted MuonStub objects
  std::ofstream refHitsFile_;       // For reference hits and restricted stubs
  
  bool filesOpened_;
  
  void openCSVFiles();
  void writeCSVHeaders();
  void exportDTPhiDigis(const boost::property_tree::ptree& procDataTree);
  void exportDTThetaDigis(const boost::property_tree::ptree& procDataTree);
  void exportCSCDigis(const boost::property_tree::ptree& procDataTree);
  void exportRPCDigis(const boost::property_tree::ptree& procDataTree);
  void exportGoldenResults(const boost::property_tree::ptree& procDataTree);
  void exportRefHitsEntry(unsigned int iProcessor,
                          unsigned int iRefHit,
                          const RefHitDef& refHitDef,
                          const std::vector<std::pair<unsigned int, MuonStubPtrs1D>>& allLayerStubs);
  void createOutputDirectory();
};

#endif /* L1T_OmtfP1_HLSDIGIEXPORTER_H_ */