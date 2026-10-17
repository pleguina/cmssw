#include "EMTFTools/ParticleGuns/interface/FlatRandomTau3MuGunProducer2.h"

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

#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

using namespace edm;
using namespace CLHEP;

FlatRandomTau3MuGunProducer2::FlatRandomTau3MuGunProducer2(const ParameterSet& pset): 
    BaseFlatGunProducer2(pset) 
{
    const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

    min_pt_Ds_ = pgun_params.getParameter<double>("MinPtDs") * GeV;
    max_pt_Ds_ = pgun_params.getParameter<double>("MaxPtDs") * GeV;
    min_invpt_Ds_ = (max_pt_Ds_ != 0.) ? 1. / max_pt_Ds_ : 1e-9;
    max_invpt_Ds_ = (min_pt_Ds_ != 0.) ? 1. / min_pt_Ds_ : 1e-9;
    random_charge_en_ = pgun_params.exists("RandomCharge") ? pgun_params.getParameter<bool>("RandomCharge") : false;

    produces<HepMCProduct>("unsmeared");
    produces<GenEventInfoProduct>();
}

FlatRandomTau3MuGunProducer2::~FlatRandomTau3MuGunProducer2() {
    // Do nothing
}

void FlatRandomTau3MuGunProducer2::produce(Event& evt, const EventSetup& es) {
    // Get particle table
    auto const& pdg_table = es.getData(pdg_table_es_token_);

    // Get random number generator
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

    if (verbosity_ > 0) {
        std::cout << " FlatRandomTau3MuGunProducer2 : Begin New Event Generation" << std::endl;
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
        // std::cout << "daughter_mass: " << daughter_mass / GeV << std::endl;

        // Random Charge
        if (random_charge_en_ && (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5)) {
            part_id_ = -part_id_;
        }

        // Get Particle Info
        double Ds_mass = 1968.35 * MeV;
        double Ds_ctau = 151.2 * um;
        double tau_mass = 1776.86  * MeV;
        double tau_ctau = 87.03 * um;

        // Calculate Ds 4-vertex, tau and nu 4-momentum
        XYZTLorentzVectorD Ds_vtx; 
        XYZTLorentzVectorD tau_4mom, nu_4mom;

        shoot_tau(
            engine, Ds_mass, Ds_ctau, tau_mass, 
            Ds_vtx, tau_4mom, nu_4mom
        );
        
        // Decay tau into 3 muons
        XYZTLorentzVectorD tau_vtx;
        XYZTLorentzVectorD mu1_4mom, mu2_4mom, mu3_4mom;

        decay_tau(
                engine, tau_mass, tau_ctau, daughter_mass, 
                tau_4mom, tau_vtx, 
                mu1_4mom, mu2_4mom, mu3_4mom
        );

        // Add Ds displacement
        tau_vtx = tau_vtx + Ds_vtx;

        // Check that 4-momentum balance out
        double tau_4p_balance = std::abs((tau_4mom - mu1_4mom - mu2_4mom - mu3_4mom).mag());

        if (tau_4p_balance > 1e-6) {
            std::cout << "tau 4-momentum not balanced: " << tau_4p_balance << std::endl;
        }

        // Convert to GeV units
        mu1_4mom = mu1_4mom / GeV;
        mu2_4mom = mu2_4mom / GeV;
        mu3_4mom = mu3_4mom / GeV;

        // Create Particles
        HepMC::FourVector mu1_p(
            mu1_4mom.Px(), mu1_4mom.Py(), mu1_4mom.Pz(), mu1_4mom.E()
        );
        HepMC::FourVector mu2_p(
            mu2_4mom.Px(), mu2_4mom.Py(), mu2_4mom.Pz(), mu2_4mom.E()
        );
        HepMC::FourVector mu3_p(
            mu3_4mom.Px(), mu3_4mom.Py(), mu3_4mom.Pz(), mu3_4mom.E()
        );

        // std::cout << "mu p1: " << mu1_4mom.Pt() << " " << mu1_4mom.Eta() << " " << mu1_4mom.Phi() << std::endl;
        // std::cout << "mu p2: " << mu2_4mom.Pt() << " " << mu2_4mom.Eta() << " " << mu2_4mom.Phi() << std::endl;
        // std::cout << "mu p3: " << mu3_4mom.Pt() << " " << mu3_4mom.Eta() << " " << mu3_4mom.Phi() << std::endl;
        // std::cout << "tau mass: " << (mu1_4mom+mu2_4mom+mu3_4mom).M() << std::endl;
        // std::cout << "tau vtx: " << tau_vtx.x() << " " << tau_vtx.y() << " " << tau_vtx.z() << std::endl;
        
        auto* mu1 = new HepMC::GenParticle(mu1_p, part_id_, 1);
        auto* mu2 = new HepMC::GenParticle(mu2_p, -part_id_, 1);
        auto* mu3 = new HepMC::GenParticle(mu3_p, -part_id_, 1);

        mu1->suggest_barcode(barcode);
        barcode++;
        mu2->suggest_barcode(barcode);
        barcode++;
        mu3->suggest_barcode(barcode);
        barcode++;
        
        // Add First Vertex
        vertex = new HepMC::GenVertex(HepMC::FourVector(
                    tau_vtx.X(), tau_vtx.Y(), tau_vtx.Z(), tau_vtx.T()
        ));
        vertex->add_particle_out(mu1);
        vertex->add_particle_out(mu2);
        vertex->add_particle_out(mu3);
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
        std::cout << " FlatRandomTau3MuGunProducer2 : Event Generation Done " << std::endl;
    }
}

void FlatRandomTau3MuGunProducer2::shoot_tau(
    CLHEP::HepRandomEngine* engine,
    const double& Ds_mass, const double& Ds_ctau, const double& tau_mass,
    XYZTLorentzVectorD& Ds_vtx, XYZTLorentzVectorD& tau_4mom, XYZTLorentzVectorD& nu_4mom
) const {
    // Calculate h 4-momentum in Lab Frame
    double Ds_eta_sign;

    if (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5) {
        Ds_eta_sign = 1.;
    } else {
        Ds_eta_sign = -1.;
    }

    double randval = CLHEP::RandFlat::shoot(engine, 0, 1);
    double Ds_pt = 1 / std::exp((1 - randval) * std::log(min_invpt_Ds_) + randval * std::log(max_invpt_Ds_));
    double Ds_eta = Ds_eta_sign * CLHEP::RandFlat::shoot(engine, min_eta_, max_eta_);
    double Ds_theta = 2 * std::atan(std::exp(-Ds_eta));
    double Ds_phi = CLHEP::RandFlat::shoot(engine, min_phi_, max_phi_);
    double Ds_p = Ds_pt / std::sin(Ds_theta);
    double Ds_e = std::hypot(Ds_p, Ds_mass);
    double Ds_px = Ds_p * std::sin(Ds_theta) * std::cos(Ds_phi);
    double Ds_py = Ds_p * std::sin(Ds_theta) * std::sin(Ds_phi);
    double Ds_pz = Ds_p * std::cos(Ds_theta);

    XYZTLorentzVectorD Ds_4mom(Ds_px, Ds_py, Ds_pz, Ds_e);

    // Decay Ds
    decay_meson(
            engine, 
            Ds_mass, Ds_ctau, tau_mass, 
            Ds_4mom, Ds_vtx, tau_4mom, nu_4mom
    );
}


void FlatRandomTau3MuGunProducer2::decay_meson(
    CLHEP::HepRandomEngine* engine,
    const double& parent_mass, const double& parent_ctau, const double& tau_mass, 
    const XYZTLorentzVectorD& parent_4mom, XYZTLorentzVectorD& decay_vtx, 
    XYZTLorentzVectorD& tau_4mom, XYZTLorentzVectorD& nu_4mom
) const {
    // Calculate parent boost vector
    auto parent_boost_beta = parent_4mom.BoostToCM();

    // Calculate daughter momentum in parent frame
    double daughter_p = (parent_mass * parent_mass - tau_mass * tau_mass) / 2 / parent_mass;
    double daughter_theta = CLHEP::RandFlat::shoot(engine, -M_PI / 2, M_PI / 2);
    double daughter_phi = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);
    double daughter_px = daughter_p * std::sin(daughter_theta) * std::cos(daughter_phi);
    double daughter_py = daughter_p * std::sin(daughter_theta) * std::sin(daughter_phi);
    double daughter_pz = daughter_p * std::cos(daughter_theta);

    tau_4mom.SetPxPyPzE(daughter_px, daughter_py, daughter_pz, std::hypot(daughter_p, tau_mass));
    nu_4mom.SetPxPyPzE(-daughter_px, -daughter_py, -daughter_pz, daughter_p);

    // Calculate decay time
    double parent_ctp = -parent_ctau * std::log(1 - CLHEP::RandFlat::shoot(engine, 0, 1));

    XYZTLorentzVectorD parent_vtx(0, 0, 0, parent_ctp); 

    // Boost 4-vectors to lab frame
    ROOT::Math::Boost boost_to_lab(-parent_boost_beta);

    parent_vtx = boost_to_lab * parent_vtx;
    tau_4mom = boost_to_lab * tau_4mom;
    nu_4mom = boost_to_lab * nu_4mom;

    // Set parent vertex as decay vertex
    decay_vtx.SetXYZT(parent_vtx.X(), parent_vtx.Y(), parent_vtx.Z(), parent_vtx.T());
}

void FlatRandomTau3MuGunProducer2::decay_tau(
    CLHEP::HepRandomEngine* engine,
    const double& parent_mass, const double& parent_ctau, const double& daughter_mass, 
    const XYZTLorentzVectorD& parent_4mom, XYZTLorentzVectorD& decay_vtx, 
    XYZTLorentzVectorD& d1_4mom, XYZTLorentzVectorD& d2_4mom, XYZTLorentzVectorD& d3_4mom
) const {
    // Calculate parent boost vector
    auto parent_boost_beta = parent_4mom.BoostToCM();

    // Calculate daughter momentums in parent frame
    double max_p = std::sqrt((parent_mass * parent_mass - 9 * daughter_mass * daughter_mass) * (parent_mass * parent_mass - daughter_mass * daughter_mass)) / parent_mass / 2;
    double max_e = std::hypot(max_p, daughter_mass);
        
    // std::cout << "max_p: " << max_p << " max_e: " << max_e << std::endl;
    
    while(true) {
        double p1 = CLHEP::RandFlat::shoot(engine, 0, max_p);
        double e1 = std::hypot(p1, daughter_mass);

        double max_e2 = std::min(parent_mass - e1 - daughter_mass, max_e);
        double min_e2 = parent_mass - e1 - max_e;
        double e2 = CLHEP::RandFlat::shoot(engine, min_e2, max_e2);
        double p2 = std::sqrt(e2 * e2 - daughter_mass * daughter_mass);

        double e3 = parent_mass - e1 - e2;
        double p3 = std::sqrt(e3 * e3 - daughter_mass * daughter_mass);

        // Calc angle p1-p2
        double arg_p1_p2 = (p1*p1 + p2*p2 - p3*p3) / 2. / p1 / p2;
        double arg_p2_p3 = (p2*p2 + p3*p3 - p1*p1) / 2. / p2 / p3;
        double arg_p3_p1 = (p3*p3 + p1*p1 - p2*p2) / 2. / p3 / p1;

        if (std::abs(arg_p1_p2) > 1) {
            std::cout << "Angle p1-p2 failed. Try again" << std::endl;
            continue;
        }

        if (std::abs(arg_p2_p3) > 1) {
            std::cout << "Angle p2-p3 failed. Try again" << std::endl;
            continue;
        }

        if (std::abs(arg_p3_p1) > 1) {
            std::cout << "Angle p3-p1 failed. Try again" << std::endl;
            continue;
        }

        double ang_p1_p2 = M_PI - std::acos(arg_p1_p2);

        // Calculate momentum vectors in orientation where p1 is along the x-axis
        ROOT::Math::XYZVectorD p1_vec(p1, 0, 0);
        ROOT::Math::XYZVectorD p2_vec(p2 * std::cos(ang_p1_p2), p2 * std::sin(ang_p1_p2), 0);
        auto p3_vec = -(p1_vec + p2_vec);

        // std::cout << "p1: " << p1_vec.Rho() << " " << p1_vec.Eta() << " " << p1_vec.Phi() << std::endl;
        // std::cout << "p2: " << p2_vec.Rho() << " " << p2_vec.Eta() << " " << p2_vec.Phi() << std::endl;
        // std::cout << "p3: " << p3_vec.Rho() << " " << p3_vec.Eta() << " " << p3_vec.Phi() << std::endl;
        
        // Calculate random rotation
        double rand_theta = CLHEP::RandFlat::shoot(engine, 0, M_PI);
        double rand_phi = CLHEP::RandFlat::shoot(engine, 0, 2 * M_PI);
        double rand_arg = CLHEP::RandFlat::shoot(engine, 0, 2 * M_PI);
        
        ROOT::Math::XYZVectorD rot_axis(
                std::sin(rand_theta) * std::cos(rand_phi), 
                std::sin(rand_theta) * std::sin(rand_phi), 
                std::cos(rand_theta) 
        );

        ROOT::Math::AxisAngle axis_angle(rot_axis, rand_arg);
        ROOT::Math::Rotation3D rand_rot(axis_angle);

        // Rotate momentum vectors
        p1_vec = rand_rot * p1_vec;
        p2_vec = rand_rot * p2_vec;
        p3_vec = rand_rot * p3_vec;

        // std::cout << "rot p1: " << p1_vec.Rho() << " " << p1_vec.Eta() << " " << p1_vec.Phi() << std::endl;
        // std::cout << "rot p2: " << p2_vec.Rho() << " " << p2_vec.Eta() << " " << p2_vec.Phi() << std::endl;
        // std::cout << "rot p3: " << p3_vec.Rho() << " " << p3_vec.Eta() << " " << p3_vec.Phi() << std::endl;
    
        d1_4mom.SetPxPyPzE(p1_vec.x(), p1_vec.y(), p1_vec.z(), e1);
        d2_4mom.SetPxPyPzE(p2_vec.x(), p2_vec.y(), p2_vec.z(), e2);
        d3_4mom.SetPxPyPzE(p3_vec.x(), p3_vec.y(), p3_vec.z(), e3);

        // exit
        break;
    }

    // Calculate decay time
    double parent_ctp = -parent_ctau * std::log(1 - CLHEP::RandFlat::shoot(engine, 0, 1));

    XYZTLorentzVectorD parent_vtx(0, 0, 0, parent_ctp); 

    // Boost 4-vectors to lab frame
    ROOT::Math::Boost boost_to_lab(-parent_boost_beta);

    parent_vtx = boost_to_lab * parent_vtx;
    d1_4mom = boost_to_lab * d1_4mom;
    d2_4mom = boost_to_lab * d2_4mom;
    d3_4mom = boost_to_lab * d3_4mom;

    // Set parent vertex as decay vertex
    decay_vtx.SetXYZT(parent_vtx.X(), parent_vtx.Y(), parent_vtx.Z(), parent_vtx.T());
}

//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomTau3MuGunProducer2);
