#include "EMTFTools/ParticleGuns/interface/FlatRandomLLPGunProducer2.h"

#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "Math/Boost.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace edm;
using namespace CLHEP;

FlatRandomLLPGunProducer2::FlatRandomLLPGunProducer2(const ParameterSet& pset): 
    BaseFlatGunProducer2(pset) 
{
    const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

    min_m_h_ = pgun_params.getParameter<double>("MinMassH") * GeV;
    max_m_h_ = pgun_params.getParameter<double>("MaxMassH") * GeV;
    min_pt_h_ = pgun_params.getParameter<double>("MinPtH") * GeV;
    max_pt_h_ = pgun_params.getParameter<double>("MaxPtH") * GeV;
    min_invpt_h_ = (max_pt_h_ != 0.) ? 1. / max_pt_h_ : 1e-9;
    max_invpt_h_ = (min_pt_h_ != 0.) ? 1. / min_pt_h_ : 1e-9;
    min_ctau_llp_ = pgun_params.getParameter<double>("MinCTauLLP") * mm;
    max_ctau_llp_ = pgun_params.getParameter<double>("MaxCTauLLP") * mm;
    llp_mass_spectrum_ = pgun_params.exists("LLPMassSpectrum") ? pgun_params.getParameter<std::string>("LLPMassSpectrum") : "flatMass";

    // Make sure min higgs mass is set to 20 GeV since min llp mass is 10 GeV 
    if (min_m_h_ < (20. * GeV)) {
        min_m_h_ = 20. * GeV;
    } 

    produces<HepMCProduct>("unsmeared");
    produces<GenEventInfoProduct>();
}

FlatRandomLLPGunProducer2::~FlatRandomLLPGunProducer2() {
    // Do nothing
}

void FlatRandomLLPGunProducer2::produce(Event& evt, const EventSetup& es) {
    // Get particle table
    auto const& pdg_table = es.getData(pdg_table_es_token_);

    // Get random number generator
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

    if (verbosity_ > 0) {
        std::cout << " FlatRandomLLPGunProducer2 : Begin New Event Generation" << std::endl;
    }

    // here re-create event (memory)
    event_ = new HepMC::GenEvent();

    // Primary vertex
    HepMC::GenVertex* vertex = nullptr;

    // Loop over particles
    int barcode = 1;

    for (unsigned int i_particle = 0; i_particle < particle_id_list_.size(); ++i_particle) {
        // Get Particle Info
        int part_id_ = particle_id_list_[i_particle];

        const HepPDT::ParticleData* particle_data = pdg_table.particle(
                HepPDT::ParticleID(abs(part_id_))
        );
        double daughter_mass = particle_data->mass().value() * GeV;

        // Calculate Higgs mass
        double max_h_mass = max_m_h_;
        double min_h_mass = min_m_h_;
        double h_mass = CLHEP::RandFlat::shoot(engine, min_h_mass, max_h_mass);

        // Calculate llp ctau
        double max_llp_ctau = max_ctau_llp_;
        double min_llp_ctau = min_ctau_llp_;
        double llp_ctau = CLHEP::RandFlat::shoot(engine, min_llp_ctau, max_llp_ctau);

        // Calculate llp mass
        double max_llp_mass = 0.5 * h_mass;
        double min_llp_mass = 10. * GeV;
        double llp_mass;

        if (llp_mass_spectrum_ == "flatInvDMass") {
            double max_dmass = h_mass - min_llp_mass;
            double min_dmass = h_mass - max_llp_mass;
            double max_invdmass = 1 / min_dmass;
            double min_invdmass = 1 / max_dmass;
            double dmass = 1 / CLHEP::RandFlat::shoot(engine, min_invdmass, max_invdmass);
            
            llp_mass = h_mass - dmass;
        } else {
            llp_mass = CLHEP::RandFlat::shoot(engine, min_llp_mass, max_llp_mass);
        }

        // Calculate Higgs 4-vertex and LLP 4-momentum
        XYZTLorentzVectorD h_vtx; 
        XYZTLorentzVectorD llp1_4mom, llp2_4mom;

        shoot_llp(
            engine, h_mass, llp_mass, 
            h_vtx, llp1_4mom, llp2_4mom
        );
        
        // Decay llp into 2 particles each
        XYZTLorentzVectorD llp1_vtx, llp2_vtx;
        XYZTLorentzVectorD llp1_out1_4mom, llp1_out2_4mom;
        XYZTLorentzVectorD llp2_out1_4mom, llp2_out2_4mom;

        decay_particle(
                engine, llp_mass, llp_ctau, daughter_mass, 
                llp1_4mom, llp1_vtx, 
                llp1_out1_4mom, llp1_out2_4mom
        );
        decay_particle(
                engine, llp_mass, llp_ctau, daughter_mass, 
                llp2_4mom, llp2_vtx, 
                llp2_out1_4mom, llp2_out2_4mom
        );

        // Add higgs displacement
        llp1_vtx = llp1_vtx + h_vtx;
        llp2_vtx = llp2_vtx + h_vtx;

        // Check that 4-momentum balance out
        double llp1_4p_balance = std::abs((llp1_4mom - llp1_out1_4mom - llp1_out2_4mom).mag());
        double llp2_4p_balance = std::abs((llp2_4mom - llp2_out1_4mom - llp2_out2_4mom).mag());

        if (llp1_4p_balance > 1e-6) {
            std::cout << "LLP1 4-momentum not balanced: " << llp1_4p_balance << std::endl;
        }

        if (llp2_4p_balance > 1e-6) {
            std::cout << "LLP2 4-momentum not balanced: " << llp2_4p_balance << std::endl;
        }

        // Convert to GeV units
        llp1_4mom = llp1_4mom / GeV;
        llp2_4mom = llp2_4mom / GeV;
        llp1_out1_4mom = llp1_out1_4mom / GeV;
        llp1_out2_4mom = llp1_out2_4mom / GeV;
        llp2_out1_4mom = llp2_out1_4mom / GeV;
        llp2_out2_4mom = llp2_out2_4mom / GeV;

        // Create Particles
        HepMC::FourVector llp1_out1_p(
            llp1_out1_4mom.Px(), llp1_out1_4mom.Py(), llp1_out1_4mom.Pz(), llp1_out1_4mom.E()
        );
        HepMC::FourVector llp1_out2_p(
            llp1_out2_4mom.Px(), llp1_out2_4mom.Py(), llp1_out2_4mom.Pz(), llp1_out2_4mom.E()
        );
        HepMC::FourVector llp2_out1_p(
            llp2_out1_4mom.Px(), llp2_out1_4mom.Py(), llp2_out1_4mom.Pz(), llp2_out1_4mom.E()
        );
        HepMC::FourVector llp2_out2_p(
            llp2_out2_4mom.Px(), llp2_out2_4mom.Py(), llp2_out2_4mom.Pz(), llp2_out2_4mom.E()
        );
        
        auto* llp1_out1 = new HepMC::GenParticle(llp1_out1_p, part_id_, 1);
        auto* llp1_out2 = new HepMC::GenParticle(llp1_out2_p, -part_id_, 1);
        auto* llp2_out1 = new HepMC::GenParticle(llp2_out1_p, part_id_, 1);
        auto* llp2_out2 = new HepMC::GenParticle(llp2_out2_p, -part_id_, 1);

        llp1_out1->suggest_barcode(barcode);
        barcode++;
        llp1_out2->suggest_barcode(barcode);
        barcode++;
        
        llp2_out1->suggest_barcode(barcode);
        barcode++;
        llp2_out2->suggest_barcode(barcode);
        barcode++;
        
        // Add First Vertex
        vertex = new HepMC::GenVertex(HepMC::FourVector(
                    llp1_vtx.X(), llp1_vtx.Y(), llp1_vtx.Z(), llp1_vtx.T()
        ));
        vertex->add_particle_out(llp1_out1);
        vertex->add_particle_out(llp1_out2);
        event_->add_vertex(vertex);

        // Add Second Vertex
        vertex = new HepMC::GenVertex(HepMC::FourVector(
                    llp2_vtx.X(), llp2_vtx.Y(), llp2_vtx.Z(), llp2_vtx.T()
        ));
        vertex->add_particle_out(llp2_out1);
        vertex->add_particle_out(llp2_out2);
        event_->add_vertex(vertex);
    }

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
        std::cout << " FlatRandomLLPGunProducer2 : Event Generation Done " << std::endl;
    }
}

void FlatRandomLLPGunProducer2::shoot_llp(
    CLHEP::HepRandomEngine* engine,
    const double& h_mass, const double& llp_mass, XYZTLorentzVectorD& h_vtx, 
    XYZTLorentzVectorD& llp1_4mom, XYZTLorentzVectorD& llp2_4mom
) const {
    // Calculate width
    double h_width = 0.027 * h_mass;
    double h_ctau = 0.19733e-15 * (GeV * m) / h_width;

    // Calculate h 4-momentum in Lab Frame
    double h_eta_sign;

    if (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5) {
        h_eta_sign = 1.;
    } else {
        h_eta_sign = -1.;
    }

    double randval = CLHEP::RandFlat::shoot(engine, 0, 1);
    double h_pt = 1 / std::exp((1 - randval) * std::log(min_invpt_h_) + randval * std::log(max_invpt_h_));
    double h_eta = h_eta_sign * CLHEP::RandFlat::shoot(engine, min_eta_, max_eta_);
    double h_theta = 2 * std::atan(std::exp(-h_eta));
    double h_phi = CLHEP::RandFlat::shoot(engine, min_phi_, max_phi_);
    double h_p = h_pt / std::sin(h_theta);
    double h_e = std::hypot(h_p, h_mass);
    double h_px = h_p * std::sin(h_theta) * std::cos(h_phi);
    double h_py = h_p * std::sin(h_theta) * std::sin(h_phi);
    double h_pz = h_p * std::cos(h_theta);

    XYZTLorentzVectorD h_4mom(h_px, h_py, h_pz, h_e);

    // Decay Higgs
    decay_particle(
            engine, 
            h_mass, h_ctau, llp_mass, h_4mom, 
            h_vtx, llp1_4mom, llp2_4mom
    );
}


void FlatRandomLLPGunProducer2::decay_particle(
    CLHEP::HepRandomEngine* engine,
    const double& parent_mass, const double& parent_ctau, const double& daughter_mass, 
    const XYZTLorentzVectorD& parent_4mom, XYZTLorentzVectorD& decay_vtx, 
    XYZTLorentzVectorD& daughter1_4mom, XYZTLorentzVectorD& daughter2_4mom
) const {
    // Calculate parent boost vector
    auto parent_boost_beta = parent_4mom.BoostToCM();

    // Calculate daughter momentum in parent frame
    double daughter_p = 0.5 * std::sqrt(parent_mass * parent_mass - 4. * daughter_mass * daughter_mass);
    double daughter_e = parent_mass / 2.;
    double daughter_theta = CLHEP::RandFlat::shoot(engine, -M_PI / 2, M_PI / 2);
    double daughter_phi = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);
    double daughter_px = daughter_p * std::sin(daughter_theta) * std::cos(daughter_phi);
    double daughter_py = daughter_p * std::sin(daughter_theta) * std::sin(daughter_phi);
    double daughter_pz = daughter_p * std::cos(daughter_theta);

    daughter1_4mom.SetPxPyPzE(daughter_px, daughter_py, daughter_pz, daughter_e);
    daughter2_4mom.SetPxPyPzE(-daughter_px, -daughter_py, -daughter_pz, daughter_e);

    // Calculate decay time
    double parent_ctp = -parent_ctau * std::log(1 - CLHEP::RandFlat::shoot(engine, 0, 1));

    XYZTLorentzVectorD parent_vtx(0, 0, 0, parent_ctp); 

    // Boost 4-vectors to lab frame
    ROOT::Math::Boost boost_to_lab(-parent_boost_beta);

    parent_vtx     = boost_to_lab * parent_vtx;
    daughter1_4mom = boost_to_lab * daughter1_4mom;
    daughter2_4mom = boost_to_lab * daughter2_4mom;

    // Set parent vertex as decay vertex
    decay_vtx.SetXYZT(parent_vtx.X(), parent_vtx.Y(), parent_vtx.Z(), parent_vtx.T());
}

//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomLLPGunProducer2);
