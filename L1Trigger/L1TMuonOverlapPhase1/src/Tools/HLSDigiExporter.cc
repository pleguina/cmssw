#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/HLSDigiExporter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iomanip>
#include <sstream>

HLSDigiExporter::HLSDigiExporter(const std::string& outputDir)
    : outputDir_(outputDir), currentEvent_(0), currentRun_(0), currentProcessor_(0), 
      currentMtfType_(l1t::omtf_pos), filesOpened_(false) {
  createOutputDirectory();
  edm::LogInfo("HLSDigiExporter") << "HLS Digi Exporter initialized. Output dir: " << outputDir_;
}

HLSDigiExporter::~HLSDigiExporter() {
  if (dtPhiDigiFile_.is_open()) dtPhiDigiFile_.close();
  if (dtThetaDigiFile_.is_open()) dtThetaDigiFile_.close();
  if (cscDigiFile_.is_open()) cscDigiFile_.close();
  if (rpcDigiFile_.is_open()) rpcDigiFile_.close();
  if (goldenResultsFile_.is_open()) goldenResultsFile_.close();
}

void HLSDigiExporter::observeEventBegin(const edm::Event& iEvent) {
  currentEvent_ = iEvent.id().event();
  currentRun_ = iEvent.id().run();
  
  if (!filesOpened_) {
    openCSVFiles();
    writeCSVHeaders();
    filesOpened_ = true;
  }
  
  edm::LogInfo("HLSDigiExporter") 
      << "Processing Event: " << currentEvent_ << " Run: " << currentRun_;
}

void HLSDigiExporter::observeEventEnd(const edm::Event& iEvent, 
                                     std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  // Flush all files
  if (dtPhiDigiFile_.is_open()) dtPhiDigiFile_.flush();
  if (dtThetaDigiFile_.is_open()) dtThetaDigiFile_.flush();
  if (cscDigiFile_.is_open()) cscDigiFile_.flush();
  if (rpcDigiFile_.is_open()) rpcDigiFile_.flush();
  if (goldenResultsFile_.is_open()) goldenResultsFile_.flush();
}

void HLSDigiExporter::observeProcesorBegin(unsigned int iProcessor, l1t::tftype mtfType) {
  currentProcessor_ = iProcessor;
  currentMtfType_ = mtfType;
}

void HLSDigiExporter::addProcesorData(std::string key, boost::property_tree::ptree& procDataTree) {
  // This is where we intercept the XML data and convert it to CSV
  if (key == "linkData") {
    // DT digis come through linkData
    exportDTPhiDigis(procDataTree);
    exportDTThetaDigis(procDataTree);
    // Also export golden results (converted stubs) that come with DT data
    exportGoldenResults(procDataTree);
  } else if (key == "CSC") {
    // CSC digis come through separate CSC key
    exportCSCDigis(procDataTree);
  } else if (key == "RPC") {
    // RPC digis come through separate RPC key
    exportRPCDigis(procDataTree);
  }
}

void HLSDigiExporter::endJob() {
  edm::LogInfo("HLSDigiExporter") << "HLS Digi Export completed. Files saved in: " << outputDir_;
}

void HLSDigiExporter::openCSVFiles() {
  std::string dtPhiFile = outputDir_ + "/dt_phi_digis.csv";
  std::string dtThetaFile = outputDir_ + "/dt_theta_digis.csv";
  std::string cscFile = outputDir_ + "/csc_digis.csv";
  std::string rpcFile = outputDir_ + "/rpc_digis.csv";
  std::string goldenFile = outputDir_ + "/golden_results.csv";
  
  dtPhiDigiFile_.open(dtPhiFile);
  dtThetaDigiFile_.open(dtThetaFile);
  cscDigiFile_.open(cscFile);
  rpcDigiFile_.open(rpcFile);
  goldenResultsFile_.open(goldenFile);
  
  if (!dtPhiDigiFile_.is_open() || !dtThetaDigiFile_.is_open() || 
      !cscDigiFile_.is_open() || !rpcDigiFile_.is_open() || !goldenResultsFile_.is_open()) {
    edm::LogError("HLSDigiExporter") << "Failed to open CSV files!";
  } else {
    edm::LogInfo("HLSDigiExporter") << "Opened CSV files for digi export";
  }
}

void HLSDigiExporter::writeCSVHeaders() {
  // DT Phi Digi header
  dtPhiDigiFile_ << "event,run,processor,tftype,dtID,whNum,scNum,stNum,slNum,quality,rpcFlag,phi,phiBend,bx\n";
  
  // DT Theta Digi header  
  dtThetaDigiFile_ << "event,run,processor,tftype,dtID,whNum,scNum,stNum,quality,rpcFlag,k,z,bx\n";
  
  // CSC Digi header - comprehensive fields (cscID already included)
  cscDigiFile_ << "event,run,processor,tftype,endcap,station,ring,chamber,layer,"
               << "trknmb,valid,quality,keywire,strip,pattern,bend,bx,mpclink,bx0,syncErr,cscID,"
               << "isRun3,quartStripBit,eighthStripBit,run3Pattern,slope,hmt,"
               << "fractionalStrip,fractionalSlope,clctPattern,stripType,bxData,type,"
               << "offset,scale,order\n";
  
  // RPC Digi header - comprehensive fields
  rpcDigiFile_ << "event,run,processor,tftype,rpcID,region,ring,station,sector,layer,subsector,roll,"
               << "strip,bx,time,coordinateX,coordinateY,deltaTime,deltaX,deltaY,"
               << "hasTime,hasX,hasY,isPseudoDigi\n";
               
  // Golden Results header - MuonStub objects after conversion
  goldenResultsFile_ << "event,run,processor,tftype,type,logicLayer,phiHw,etaHw,qualityHw,"
                     << "phiBHw,bx,timing,r,detId\n";
}

void HLSDigiExporter::exportDTPhiDigis(const boost::property_tree::ptree& procDataTree) {
  try {
    // Iterate through all dtP2PhiDigi entries in the procDataTree
    for (const auto& child : procDataTree) {
      if (child.first == "dtP2PhiDigi" || child.first == "dtPhiDigi") {
        const auto& digi = child.second;
        
        dtPhiDigiFile_ << currentEvent_ << ","
                       << currentRun_ << ","
                       << currentProcessor_ << ","
                       << (currentMtfType_ == l1t::omtf_neg ? "NEG" : 
                          (currentMtfType_ == l1t::omtf_pos ? "POS" : "BARREL")) << ","
                       << digi.get<unsigned int>("<xmlattr>.dtID", 0) << ","
                       << digi.get<int>("<xmlattr>.whNum") << ","
                       << digi.get<int>("<xmlattr>.scNum") << ","
                       << digi.get<int>("<xmlattr>.stNum") << ","
                       << digi.get<int>("<xmlattr>.slNum") << ","
                       << digi.get<int>("<xmlattr>.quality") << ","
                       << digi.get<int>("<xmlattr>.rpcFlag") << ","
                       << digi.get<int>("<xmlattr>.phi") << ","
                       << digi.get<int>("<xmlattr>.phiBend") << ","
                       << digi.get<int>("<xmlattr>.bx", 0) << "\n"; // bx might not always be present
      }
    }
  } catch (const std::exception& e) {
    edm::LogWarning("HLSDigiExporter") << "Error exporting DT Phi digis: " << e.what();
  }
}

void HLSDigiExporter::exportDTThetaDigis(const boost::property_tree::ptree& procDataTree) {
  try {
    for (const auto& child : procDataTree) {
      if (child.first == "dtP2ThDigi" || child.first == "dtThDigi") {
        const auto& digi = child.second;
        
        dtThetaDigiFile_ << currentEvent_ << ","
                         << currentRun_ << ","
                         << currentProcessor_ << ","
                         << (currentMtfType_ == l1t::omtf_neg ? "NEG" : 
                            (currentMtfType_ == l1t::omtf_pos ? "POS" : "BARREL")) << ","
                         << digi.get<unsigned int>("<xmlattr>.dtID", 0) << ","
                         << digi.get<int>("<xmlattr>.whNum") << ","
                         << digi.get<int>("<xmlattr>.scNum") << ","
                         << digi.get<int>("<xmlattr>.stNum") << ","
                         << digi.get<int>("<xmlattr>.quality") << ","
                         << digi.get<int>("<xmlattr>.rpcFlag") << ","
                         << digi.get<int>("<xmlattr>.k") << ","
                         << digi.get<int>("<xmlattr>.z") << ","
                         << digi.get<int>("<xmlattr>.bx", 0) << "\n";
      }
    }
  } catch (const std::exception& e) {
    edm::LogWarning("HLSDigiExporter") << "Error exporting DT Theta digis: " << e.what();
  }
}

void HLSDigiExporter::exportCSCDigis(const boost::property_tree::ptree& procDataTree) {
  try {
    for (const auto& child : procDataTree) {
      if (child.first == "cscDigi" || child.first == "cscCorrelatedLCTDigi") {
        const auto& digi = child.second;
        
        cscDigiFile_ << currentEvent_ << ","
                     << currentRun_ << ","
                     << currentProcessor_ << ","
                     << (currentMtfType_ == l1t::omtf_neg ? "NEG" : 
                        (currentMtfType_ == l1t::omtf_pos ? "POS" : "BARREL")) << ","
                     // Detector ID fields
                     << digi.get<int>("<xmlattr>.endcap", 0) << ","
                     << digi.get<int>("<xmlattr>.station", 0) << ","
                     << digi.get<int>("<xmlattr>.ring", 0) << ","
                     << digi.get<int>("<xmlattr>.chamber", 0) << ","
                     << digi.get<int>("<xmlattr>.layer", 0) << ","
                     // Core LCT data fields
                     << digi.get<int>("<xmlattr>.trknmb", 0) << ","
                     << digi.get<bool>("<xmlattr>.valid", false) << ","
                     << digi.get<int>("<xmlattr>.quality", 0) << ","
                     << digi.get<int>("<xmlattr>.keywire", 0) << ","
                     << digi.get<int>("<xmlattr>.strip", 0) << ","
                     << digi.get<int>("<xmlattr>.pattern", 0) << ","
                     << digi.get<int>("<xmlattr>.bend", 0) << ","
                     << digi.get<int>("<xmlattr>.bx", 0) << ","
                     << digi.get<int>("<xmlattr>.mpclink", 0) << ","
                     << digi.get<int>("<xmlattr>.bx0", 0) << ","
                     << digi.get<int>("<xmlattr>.syncErr", 0) << ","
                     << digi.get<int>("<xmlattr>.cscID", 0) << ","
                     // Run-3 specific fields
                     << digi.get<bool>("<xmlattr>.isRun3", false) << ","
                     << digi.get<bool>("<xmlattr>.quartStripBit", false) << ","
                     << digi.get<bool>("<xmlattr>.eighthStripBit", false) << ","
                     << digi.get<int>("<xmlattr>.run3Pattern", 0) << ","
                     << digi.get<int>("<xmlattr>.slope", 0) << ","
                     << digi.get<int>("<xmlattr>.hmt", 0) << ","
                     // Additional computed fields
                     << digi.get<float>("<xmlattr>.fractionalStrip", 0.0) << ","
                     << digi.get<float>("<xmlattr>.fractionalSlope", 0.0) << ","
                     << digi.get<int>("<xmlattr>.clctPattern", 0) << ","
                     << digi.get<int>("<xmlattr>.stripType", 0) << ","
                     << digi.get<int>("<xmlattr>.bxData", 0) << ","
                     << digi.get<int>("<xmlattr>.type", 0) << ","
                     // CSC conversion parameters
                     << digi.get<int>("<xmlattr>.offset", 0) << ","
                     << digi.get<double>("<xmlattr>.scale", 0.0) << ","
                     << digi.get<int>("<xmlattr>.order", 0) << "\n";
      }
    }
  } catch (const std::exception& e) {
    edm::LogWarning("HLSDigiExporter") << "Error exporting CSC digis: " << e.what();
  }
}

void HLSDigiExporter::exportRPCDigis(const boost::property_tree::ptree& procDataTree) {
  try {
    for (const auto& child : procDataTree) {
      if (child.first == "rpcDigi") {
        const auto& digi = child.second;
        
        rpcDigiFile_ << currentEvent_ << ","
                     << currentRun_ << ","
                     << currentProcessor_ << ","
                     << (currentMtfType_ == l1t::omtf_neg ? "NEG" : 
                        (currentMtfType_ == l1t::omtf_pos ? "POS" : "BARREL")) << ","
                     // Detector ID 
                     << digi.get<unsigned int>("<xmlattr>.rpcID", 0) << ","
                     // Detector ID fields
                     << digi.get<int>("<xmlattr>.region", 0) << ","
                     << digi.get<int>("<xmlattr>.ring", 0) << ","
                     << digi.get<int>("<xmlattr>.station", 0) << ","
                     << digi.get<int>("<xmlattr>.sector", 0) << ","
                     << digi.get<int>("<xmlattr>.layer", 0) << ","
                     << digi.get<int>("<xmlattr>.subsector", 0) << ","
                     << digi.get<int>("<xmlattr>.roll", 0) << ","
                     // Core RPC digi fields
                     << digi.get<int>("<xmlattr>.strip", 0) << ","
                     << digi.get<int>("<xmlattr>.bx", 0) << ","
                     // Timing and spatial information
                     << digi.get<double>("<xmlattr>.time", 0.0) << ","
                     << digi.get<double>("<xmlattr>.coordinateX", 0.0) << ","
                     << digi.get<double>("<xmlattr>.coordinateY", 0.0) << ","
                     << digi.get<double>("<xmlattr>.deltaTime", 0.0) << ","
                     << digi.get<double>("<xmlattr>.deltaX", 0.0) << ","
                     << digi.get<double>("<xmlattr>.deltaY", 0.0) << ","
                     // Flags for available information
                     << digi.get<bool>("<xmlattr>.hasTime", false) << ","
                     << digi.get<bool>("<xmlattr>.hasX", false) << ","
                     << digi.get<bool>("<xmlattr>.hasY", false) << ","
                     << digi.get<bool>("<xmlattr>.isPseudoDigi", false) << "\n";
      }
    }
  } catch (const std::exception& e) {
    edm::LogWarning("HLSDigiExporter") << "Error exporting RPC digis: " << e.what();
  }
}

void HLSDigiExporter::exportGoldenResults(const boost::property_tree::ptree& procDataTree) {
  try {
    for (const auto& child : procDataTree) {
      if (child.first == "goldenStub") {
        const auto& stub = child.second;
        
        goldenResultsFile_ << currentEvent_ << ","
                           << currentRun_ << ","
                           << currentProcessor_ << ","
                           << (currentMtfType_ == l1t::omtf_neg ? "NEG" : 
                              (currentMtfType_ == l1t::omtf_pos ? "POS" : "BARREL")) << ","
                           << stub.get<int>("<xmlattr>.type", 0) << ","
                           << stub.get<int>("<xmlattr>.logicLayer", 0) << ","
                           << stub.get<int>("<xmlattr>.phiHw", 0) << ","
                           << stub.get<int>("<xmlattr>.etaHw", 0) << ","
                           << stub.get<int>("<xmlattr>.qualityHw", 0) << ","
                           << stub.get<int>("<xmlattr>.phiBHw", 0) << ","
                           << stub.get<int>("<xmlattr>.bx", 0) << ","
                           << stub.get<int>("<xmlattr>.timing", 0) << ","
                           << stub.get<int>("<xmlattr>.r", 0) << ","
                           << stub.get<unsigned long>("<xmlattr>.detId", 0) << "\n";
      }
    }
  } catch (const std::exception& e) {
    edm::LogWarning("HLSDigiExporter") << "Error exporting golden results: " << e.what();
  }
}

void HLSDigiExporter::createOutputDirectory() {
  boost::filesystem::path dir(outputDir_);
  if (!boost::filesystem::exists(dir)) {
    boost::filesystem::create_directories(dir);
    edm::LogInfo("HLSDigiExporter") << "Created output directory: " << outputDir_;
  }
}