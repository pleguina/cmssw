// RpcBarrelKeyToRDumper
//
// Barrel counterpart to RpcKeyToRDumper.cc. Dumps eta/phi/r for every
// real barrel RPC roll+strip (RB1in/RB1out/RB2in/RB2out/RB3), using the
// same OmtfAngleConverter calls the production addRPCstub() path uses.
//
// r is actually a hardcoded per-layer constant in OMTFinputMaker.cc
// (RB1in=413.675, RB1out=448.675, RB2in=494.975, RB2out=529.975,
// RB3=602.150 mm), so the r dumped here is only a cross-check against
// those constants, not the value the firmware interface should use.
//
// Not filtered to any wheel/sector subset -- dumps all rolls, same as
// RpcKeyToRDumper.cc.

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfAngleConverter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/ProcConfigurationBase.h"

#include <cmath>
#include <fstream>
#include <vector>

// Same stub as RpcKeyToRDumper.cc's RpcMinimalProcConfig.
class RpcBarrelMinimalProcConfig : public ProcConfigurationBase {
public:
  unsigned int nPhiBins() const override { return 5400; }
  unsigned int nProcessors() const override { return 3; }
  double hwPtToGev(int) const override { return 0; }
  int ptGevToHw(double) const override { return 0; }
  int getProcScalePhi(double, double) const override { return 0; }
  int etaToHwEta(double eta) const override { return (int)std::lround(std::abs(eta) / etaUnit_); }
  int mb1W2Eta() const override { return 0; }
  int mb2W2Eta() const override { return 0; }
  int mb3W2Eta() const override { return 0; }
  int mb4W2Eta() const override { return 0; }
  unsigned int nLayers() const override { return 18; }
  bool isBendingLayer(unsigned int) const override { return false; }
};

class RpcBarrelKeyToRDumper : public edm::one::EDAnalyzer<> {
public:
  explicit RpcBarrelKeyToRDumper(const edm::ParameterSet& cfg)
      : outPath_(cfg.getParameter<std::string>("outPath")),
        cscGeometryEsToken_(esConsumes<CSCGeometry, MuonGeometryRecord>()),
        rpcGeometryEsToken_(esConsumes<RPCGeometry, MuonGeometryRecord>()),
        dtGeometryEsToken_(esConsumes<DTGeometry, MuonGeometryRecord>()) {}

  void analyze(const edm::Event&, const edm::EventSetup& es) override {
    const RPCGeometry& rpcGeometry = es.getData(rpcGeometryEsToken_);

    RpcBarrelMinimalProcConfig config;
    config.setStubEtaEncoding(ProcConfigurationBase::StubEtaEncoding::valueP1Scale);
    MuonGeometryTokens tokens{rpcGeometryEsToken_, cscGeometryEsToken_, dtGeometryEsToken_};
    OmtfAngleConverter angleConverter;
    angleConverter.checkAndUpdateGeometry(es, &config, tokens);

    std::ofstream out(outPath_);
    out << "wheel,sector,subsector,station,layer,roll,strip,etaHw,r_check,phiBin\n";

    for (const RPCRoll* roll : rpcGeometry.rolls()) {
      if (!roll)
        continue;
      RPCDetId id = roll->id();

      // Barrel only, stations 1-3. Station 3 (RB3) has no layer 2.
      if (id.region() != 0)
        continue;
      if (id.station() < 1 || id.station() > 3)
        continue;
      if (id.station() == 3 && id.layer() != 1)
        continue;

      const int nstrips = roll->nstrips();
      for (int strip = 1; strip <= nstrips; ++strip) {
        float r = 0;
        int etaHw = angleConverter.getGlobalEtaRpc(id.rawId(), (unsigned int)strip, r);
        // phiZero=0: raw angle bin, same convention as RpcKeyToRDumper.cc.
        int phiBin = angleConverter.getProcessorPhi(
            /*phiZero=*/0, l1t::tftype::omtf_pos, id, (unsigned int)strip, (unsigned int)strip);

        out << id.ring() << "," << id.sector() << "," << id.subsector() << "," << id.station() << "," << id.layer()
            << "," << id.roll() << "," << strip << "," << etaHw << "," << r << "," << phiBin << "\n";
      }
    }
    out.close();
    edm::LogInfo("RpcBarrelKeyToRDumper") << "Wrote " << outPath_;
  }

private:
  std::string outPath_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryEsToken_;
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeometryEsToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryEsToken_;
};

DEFINE_FWK_MODULE(RpcBarrelKeyToRDumper);
