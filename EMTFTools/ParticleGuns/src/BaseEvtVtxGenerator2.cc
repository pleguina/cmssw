#include "EMTFTools/ParticleGuns/interface/BaseEvtVtxGenerator2.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//#include "HepMC/GenEvent.h"
//#include "CLHEP/Vector/ThreeVector.h"
//#include "HepMC/SimpleVector.h"

using namespace edm;
using namespace CLHEP;

BaseEvtVtxGenerator2::BaseEvtVtxGenerator2(const ParameterSet& pset) {
    Service<RandomNumberGenerator> rng;

    if (!rng.isAvailable()) {
        throw cms::Exception("Configuration") << "The BaseEvtVtxGenerator2 requires the RandomNumberGeneratorService\n"
            "which is not present in the configuration file. \n"
            "You must add the service\n"
            "in the configuration file or remove the modules that require it.";
    }

    source_token_ = consumes<edm::HepMCProduct>(pset.getParameter<edm::InputTag>("src"));
    produces<edm::HepMCProduct>();
}

BaseEvtVtxGenerator2::~BaseEvtVtxGenerator2() {
    // Do nothing
}

void BaseEvtVtxGenerator2::produce(Event& evt, const EventSetup&) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

    Handle<HepMCProduct> HepUnsmearedMCEvt;

    evt.getByToken(source_token_, HepUnsmearedMCEvt);

    // Copy the HepMC::GenEvent
    HepMC::GenEvent* genevt = new HepMC::GenEvent(*HepUnsmearedMCEvt->GetEvent());
    std::unique_ptr<edm::HepMCProduct> HepMCEvt(new edm::HepMCProduct(genevt));

    // generate new vertex & apply the shift
    HepMCEvt->applyVtxGen(newVertex(engine));

    //HepMCEvt->LorentzBoost( 0., 142.e-6 );
    HepMCEvt->boostToLab(GetInvLorentzBoost(), "vertex");
    HepMCEvt->boostToLab(GetInvLorentzBoost(), "momentum");

    evt.put(std::move(HepMCEvt));

    return;
}
