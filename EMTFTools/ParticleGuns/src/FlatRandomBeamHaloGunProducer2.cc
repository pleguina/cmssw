#include "EMTFTools/ParticleGuns/interface/FlatRandomBeamHaloGunProducer2.h"

#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace edm;
using namespace CLHEP;

FlatRandomBeamHaloGunProducer2::FlatRandomBeamHaloGunProducer2(const ParameterSet& pset): 
    BaseFlatGunProducer2(pset) 
{
    const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

    min_p_ = pgun_params.getParameter<double>("MinP");
    max_p_ = pgun_params.getParameter<double>("MaxP");
    min_vr_ = pgun_params.getParameter<double>("MinVR") * cm;
    max_vr_ = pgun_params.getParameter<double>("MaxVR") * cm;
    min_vz_ = pgun_params.getParameter<double>("MinVZ") * cm;
    max_vz_ = pgun_params.getParameter<double>("MaxVZ") * cm;
    min_tz_ = pgun_params.getParameter<double>("MinTZ") * cm;
    max_tz_ = pgun_params.getParameter<double>("MaxTZ") * cm;
    random_charge_en_ = pgun_params.exists("RandomCharge") ? pgun_params.getParameter<bool>("RandomCharge") : false;
    p_spectrum_ = pgun_params.exists("PSpectrum") ? pgun_params.getParameter<std::string>("PSpectrum") : "flatP";

    produces<HepMCProduct>("unsmeared");
    produces<GenEventInfoProduct>();
}

FlatRandomBeamHaloGunProducer2::~FlatRandomBeamHaloGunProducer2() {
    // Do nothing
}

void FlatRandomBeamHaloGunProducer2::produce(Event& evt, const EventSetup& es) {
    // Get particle table
    auto const& pdg_table = es.getData(pdg_table_es_token_);

    // Get random number generator
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

    if (verbosity_ > 0) {
        std::cout << " FlatRandomBeamHaloGunProducer2 : Begin New Event Generation" << std::endl;
    }

    // here re-create event (memory)
    event_ = new HepMC::GenEvent();

    // Primary vertex
    HepMC::GenVertex* vertex = nullptr;

    // Loop over particles
    int barcode = 1;

    for (unsigned int i_particle = 0; i_particle < particle_id_list_.size(); ++i_particle) {
        // Random mom
        double mom = 0.;

        {
            double rand_val = CLHEP::RandFlat::shoot(engine, 0., 1.);

            if (p_spectrum_ == "flatP") {
                mom = min_p_ + rand_val * (max_p_ - min_p_);
            }
        }

        // Particle ID
        int part_id_ = particle_id_list_[i_particle];

        if (random_charge_en_ && (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5)) {
            part_id_ = -part_id_;
        }
    
        // Random detector side
        double z_sign;

        if (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5) {
            z_sign = 1.;
        } else {
            z_sign = -1.;
        }

        // Calculate vertex
        double vr = CLHEP::RandFlat::shoot(engine, min_vr_, max_vr_);
        double vrot = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);
        double vx = vr * std::sin(vrot);
        double vy = vr * std::cos(vrot);
        double vz = z_sign * CLHEP::RandFlat::shoot(engine, min_vz_, max_vz_);
        double vt = 0;

        vertex = new HepMC::GenVertex(HepMC::FourVector(vx, vy, vz, vt));

        // Calculate target
        double tr = CLHEP::RandGaussQ::shoot(engine, vr, 1. * cm);
        double tx = tr * std::sin(vrot);
        double ty = tr * std::cos(vrot);
        double tz = z_sign * CLHEP::RandFlat::shoot(engine, min_tz_, max_tz_);

        // Calculate momentum unit vector
        double dx = tx - vx;
        double dy = ty - vy;
        double dz = tz - vz;
        double magd = std::hypot(dx, dy, dz);
        
        dx = dx / magd;
        dy = dy / magd;
        dz = dz / magd;

        double dr = std::hypot(dx, dy);
        double theta = std::atan2(dr, dz);
        double phi = std::atan2(dy, dx);

        // Calculate data
        const HepPDT::ParticleData* particle_data = pdg_table.particle(
                HepPDT::ParticleID(abs(part_id_))
        );
        double mass = particle_data->mass().value();
        double pt = mom * std::sin(theta);
        double px = pt * std::cos(phi);
        double py = pt * std::sin(phi);
        double pz = mom * std::cos(theta);
        double energy2 = mom * mom + mass * mass;
        double energy = std::sqrt(energy2);

        HepMC::FourVector p(px, py, pz, energy);
        HepMC::GenParticle* particle = new HepMC::GenParticle(p, part_id_, 1);
        
        particle->suggest_barcode(barcode);
        barcode++;

        vertex->add_particle_out(particle);

        if (add_anti_particle_en_) {
            HepMC::FourVector anti_p(-px, -py, -pz, energy);
            int anti_part_id_ = -part_id_;

            if (part_id_ == 22 || part_id_ == 23) {
                anti_part_id_ = part_id_;
            }

            HepMC::GenParticle* anti_particle = new HepMC::GenParticle(anti_p, anti_part_id_, 1);
            anti_particle->suggest_barcode(barcode);
            barcode++;

            vertex->add_particle_out(anti_particle);
        }
    }

    event_->add_vertex(vertex);
    event_->set_event_number(evt.id().event());
    event_->set_signal_process_id(20);

    if (verbosity_ > 0) {
        event_->print();
    }

    std::unique_ptr<HepMCProduct> BProduct(new HepMCProduct(event_));
    evt.put(std::move(BProduct), "unsmeared");

    std::unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(event_));
    evt.put(std::move(genEventInfo));

    if (verbosity_ > 0) {
        std::cout << " FlatRandomBeamHaloGunProducer2 : Event Generation Done " << std::endl;
    }
}

//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomBeamHaloGunProducer2);
