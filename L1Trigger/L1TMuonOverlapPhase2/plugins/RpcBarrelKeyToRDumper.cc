// RpcBarrelKeyToRDumper
//
// Standalone geometry dumper for WP13 barrel follow-up (RPC barrel
// eta/phi, omtf-firmware). Sibling of RpcKeyToRDumper.cc (endcap-only,
// hard-filters region()==0 out) -- this one exhaustively iterates every
// REAL RPC barrel roll (region()==0) in the three stations OMTF's real
// logic-layer scheme reads (RB1in/RB1out/RB2in/RB2out/RB3 -- station 1-3,
// layer 1/2 for stations 1-2, layer 1 only for station 3, per
// hwToLogicLayer_0x0210.xml logic layers 10-14) and every valid strip in
// each roll, calling the exact same production code path
// RpcDigiToStubsConverterOmtf::addRPCstub() (OMTFinputMaker.cc) uses when
// building a real barrel RPC stub:
//   - eta: OmtfAngleConverter::getGlobalEtaRpc(rawid, strip, r) -- the SAME
//     function used for endcap, generic over region()==0 (confirmed by
//     reading OmtfAngleConverter.cc directly: roll->toGlobal/centreOfStrip
//     are region-agnostic; only `r` itself is region-conditional).
//   - phi: OmtfAngleConverter::getProcessorPhi(phiZero=0, part, rollId,
//     strip, strip) -- also the SAME RPCDetId-overload function endcap
//     uses, equally region-agnostic.
//   - r: OMTFinputMaker.cc's own addRPCstub() OVERWRITES whatever
//     getGlobalEtaRpc returned for r with a HARDCODED per-logic-layer
//     constant for barrel (RB1in=413.675, RB1out=448.675, RB2in=494.975,
//     RB2out=529.975, RB3=602.150 mm) -- confirmed by direct read of
//     OMTFinputMaker.cc:214-223. So r needs NO geometry dump at all for
//     barrel; it is a compile-time constant already known from source.
//     Dumped here anyway (gp.perp(), real per-roll geometry) purely as an
//     independent cross-check against those hardcoded constants, NOT as
//     the value the eventual barrel interface module should use.
//
// Dumps (wheel,sector,subsector,station,layer,roll,strip,etaHw,r_check,
// phiBin) to CSV. wheel is NOT filtered to OMTF's real-connected wheels
// (ring==+-2 only, confirmed in OMTFinputMaker.cc's acceptDigi()) --
// dumped exhaustively across all 5 real wheels, same "exhaustive, not
// pre-filtered to a downstream restriction" philosophy RpcKeyToRDumper.cc
// itself already uses (it doesn't pre-filter by sector/subsector either).

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

// Identical stub to RpcKeyToRDumper.cc's RpcMinimalProcConfig -- see that
// file's own header comment for why each pure virtual is stubbed the way
// it is. Same real etaToHwEta() implementation (valueP1Scale formula).
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

      // hwToLogicLayer_0x0210.xml maps logic layers 10-14 to RB1in/RB1out/
      // RB2in/RB2out/RB3 -- barrel (region==0), station 1-3. Station 1/2
      // have two real layers (1=in, 2=out); station 3 has only layer 1
      // (RB3 -- OMTFinputMaker.cc's own getInputNumber() forces iRoll=1,
      // nInputsPerSector=2 for station==3, confirmed by direct read).
      if (id.region() != 0)
        continue;
      if (id.station() < 1 || id.station() > 3)
        continue;
      if (id.station() == 3 && id.layer() != 1)
        continue;  // RB3 has no real "layer 2"; skip to avoid a duplicate/non-existent entry

      const int nstrips = roll->nstrips();
      for (int strip = 1; strip <= nstrips; ++strip) {
        float r = 0;
        int etaHw = angleConverter.getGlobalEtaRpc(id.rawId(), (unsigned int)strip, r);
        // phiZero=0: raw, processor-independent angle bin, same convention
        // RpcKeyToRDumper.cc uses (real per-processor phiZero is a cheap
        // closed-form subtraction applied in firmware, not baked in here).
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
