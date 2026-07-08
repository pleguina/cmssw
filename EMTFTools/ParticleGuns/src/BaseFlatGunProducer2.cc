/*
 *  \author Julia Yarba
 */

#include "EMTFTools/ParticleGuns/interface/BaseFlatGunProducer2.h"

#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

using namespace edm;
using namespace CLHEP;

BaseFlatGunProducer2::BaseFlatGunProducer2(const ParameterSet& pset): 
    event_(nullptr),
    pdg_table_es_token_(esConsumes<HepPDT::ParticleDataTable, PDTRecord>())
{
    Service<RandomNumberGenerator> rng;

    if (!rng.isAvailable()) {
        throw cms::Exception("Configuration")
            << "The RandomNumberProducer module requires the RandomNumberGeneratorService\n"
            "which appears to be absent.  Please add that service to your configuration\n"
            "or remove the modules that require it.";
    }

    const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

    // although there's the method ParameterSet::empty(),
    // it looks like it's NOT even necessary to check if it is,
    // before trying to extract parameters - if it is empty,
    // the default values seem to be taken
    particle_id_list_ = pgun_params.getParameter<std::vector<int> >("PartID");
    min_eta_ = pgun_params.getParameter<double>("MinEta");
    max_eta_ = pgun_params.getParameter<double>("MaxEta");
    min_phi_ = pgun_params.getParameter<double>("MinPhi");
    max_phi_ = pgun_params.getParameter<double>("MaxPhi");
    verbosity_ = pset.getUntrackedParameter<int>("Verbosity", 0);
    add_anti_particle_en_ = pset.getParameter<bool>("AddAntiParticle");
    produces<GenRunInfoProduct, Transition::EndRun>();
}

BaseFlatGunProducer2::~BaseFlatGunProducer2() {
    // Do Nothing
}

void BaseFlatGunProducer2::beginRun(const Run& run, const EventSetup& es) {
    // Do Nothing
}

void BaseFlatGunProducer2::endRun(const Run& run, const EventSetup& es) {
    // Do Nothing
}

void BaseFlatGunProducer2::endRunProduce(Run& run, const EventSetup& es) {
    // just create an empty product
    // to keep the EventContent definitions happy
    // later on we might put the info into the run info that this is a PGun
    std::unique_ptr<GenRunInfoProduct> genRunInfo(new GenRunInfoProduct());
    run.put(std::move(genRunInfo));
}
