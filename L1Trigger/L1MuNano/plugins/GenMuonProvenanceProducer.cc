/*
 * GenMuonProvenanceProducer.cc
 *
 * Produces provenance ValueMaps for generator particles, keyed to the input
 * genParticles collection:
 *   - parentPdgId      : immediate mother PDG id (0 if unavailable)
 *   - isDecayInFlight  : 1 if first non-muon ancestor is pi/K (charged/neutral)
 *                        0 otherwise
 *
 * These maps are attached to the GenMuon Nano table via ExtVar.
 */

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>
#include <cmath>

namespace {
bool isDifParentPdgId(const int pdgIdAbs) {
  return pdgIdAbs == 211 || pdgIdAbs == 321 || pdgIdAbs == 130 || pdgIdAbs == 310;
}
}  // namespace

class GenMuonProvenanceProducer : public edm::global::EDProducer<> {
public:
  explicit GenMuonProvenanceProducer(const edm::ParameterSet& cfg)
      : srcToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("src"))) {
    produces<edm::ValueMap<int>>("parentPdgId");
    produces<edm::ValueMap<int>>("isDecayInFlight");
  }

  ~GenMuonProvenanceProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("genParticles"));
    descriptions.add("GenMuonProvenanceProducer", desc);
  }

private:
  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto genParticles = event.getHandle(srcToken_);
    const size_t n = genParticles->size();

    std::vector<int> parentPdgId(n, 0);
    std::vector<int> isDecayInFlight(n, 0);

    for (size_t i = 0; i < n; ++i) {
      const reco::GenParticle& gp = (*genParticles)[i];

      const reco::Candidate* mother = (gp.numberOfMothers() > 0) ? gp.mother(0) : nullptr;
      parentPdgId[i] = mother ? mother->pdgId() : 0;

      const reco::Candidate* ancestor = mother;
      int guard = 0;
      while (ancestor && std::abs(ancestor->pdgId()) == 13 && guard < 20) {
        ancestor = (ancestor->numberOfMothers() > 0) ? ancestor->mother(0) : nullptr;
        ++guard;
      }

      if (ancestor) {
        isDecayInFlight[i] = isDifParentPdgId(std::abs(ancestor->pdgId())) ? 1 : 0;
      }
    }

    auto putMap = [&](std::vector<int>& vals, const std::string& label) {
      auto vm = std::make_unique<edm::ValueMap<int>>();
      edm::ValueMap<int>::Filler filler(*vm);
      filler.insert(genParticles, vals.begin(), vals.end());
      filler.fill();
      event.put(std::move(vm), label);
    };

    putMap(parentPdgId, "parentPdgId");
    putMap(isDecayInFlight, "isDecayInFlight");
  }

  const edm::EDGetTokenT<reco::GenParticleCollection> srcToken_;
};

DEFINE_FWK_MODULE(GenMuonProvenanceProducer);
