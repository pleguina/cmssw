#include "EMTFTools/ParticleGuns/interface/FlatRandomPtGunProducer2.h"

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

FlatRandomPtGunProducer2::FlatRandomPtGunProducer2(const ParameterSet& pset): 
    BaseFlatGunProducer2(pset) 
{
    const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

    min_pt_ = pgun_params.getParameter<double>("MinPt");
    max_pt_ = pgun_params.getParameter<double>("MaxPt");
    min_dxy_ = pgun_params.getParameter<double>("MinDxy") * cm;
    max_dxy_ = pgun_params.getParameter<double>("MaxDxy") * cm;
    min_one_over_pt_ = (max_pt_ != 0.) ? 1. / max_pt_ : 1e-9;
    max_one_over_pt_ = (min_pt_ != 0.) ? 1. / min_pt_ : 1e-9;
    min_one_over_dxy_ = (max_dxy_ != 0.) ? 1. / max_dxy_ : 1e-9;
    max_one_over_dxy_ = (min_dxy_ != 0.) ? 1. / min_dxy_ : 1e-9;

    max_lxy_ = (pgun_params.exists("MaxLxy") ? pgun_params.getParameter<double>("MaxLxy") : 300.) * cm;
    random_charge_en_ = pgun_params.exists("RandomCharge") ? pgun_params.getParameter<bool>("RandomCharge") : false;
    pt_spectrum_ = pgun_params.exists("PtSpectrum") ? pgun_params.getParameter<std::string>("PtSpectrum") : "flatPt";
    vertex_spectrum_ = pgun_params.exists("VertexSpectrum") ? pgun_params.getParameter<std::string>("VertexSpectrum") : "none";

    produces<HepMCProduct>("unsmeared");
    produces<GenEventInfoProduct>();
}

FlatRandomPtGunProducer2::~FlatRandomPtGunProducer2() {
    // Do nothing
}

void FlatRandomPtGunProducer2::produce(Event& evt, const EventSetup& es) {
    // Get particle table
    auto const& pdg_table = es.getData(pdg_table_es_token_);

    // Get random number generator
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(evt.streamID());

    if (verbosity_ > 0) {
        std::cout << " FlatRandomPtGunProducer2 : Begin New Event Generation" << std::endl;
    }

    // here re-create event (memory)
    event_ = new HepMC::GenEvent();

    // Primary vertex
    HepMC::GenVertex* vertex = nullptr;

    // Loop over particles
    int barcode = 1;

    for (unsigned int i_particle = 0; i_particle < particle_id_list_.size(); ++i_particle) {
        // Random pT
        double pt = 0.;

        {
            double rand_val = CLHEP::RandFlat::shoot(engine, 0., 1.);

            if (pt_spectrum_ == "flatPt") {
                pt = min_pt_ + rand_val * (max_pt_ - min_pt_);
            } else if (pt_spectrum_ == "flatOneOverPt") {
                pt = min_one_over_pt_ + rand_val * (max_one_over_pt_ - min_one_over_pt_);

                if (pt != 0.) {
                    pt = 1. / pt;           
                }
            } else if (pt_spectrum_ == "flatOneOverPtCMS") {
                pt = std::exp((1. - rand_val) * std::log(min_one_over_pt_) + rand_val * std::log(max_one_over_pt_));

                if (pt != 0.){
                    pt = 1. / pt;
                }
            }
        }

        // Particle ID
        int part_id_ = particle_id_list_[i_particle];

        if (random_charge_en_ && (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5)) {
            part_id_ = -part_id_;
        }
    
        // Random Eta/Phi
        double eta = CLHEP::RandFlat::shoot(engine, min_eta_, max_eta_);
        double phi = CLHEP::RandFlat::shoot(engine, min_phi_, max_phi_);

        // Random Vertex
        int charge = (part_id_ > 0) ? -1 : 1;

        if (vertex == nullptr) {
            if (vertex_spectrum_ == "flatVtx") {
                // Try to find vertex until eta_st2 is found
                double eta_st2 = eta;
                eta = 0;

                for (int itry = 0; (itry < 100) && (eta == 0); ++itry) {
                    // Delete Vertex
                    if (vertex != nullptr) {
                        delete vertex; 
                    }

                    // Calculate vertex
                    double lxy = CLHEP::RandFlat::shoot(engine, 0., max_lxy_);
                    double rot = CLHEP::RandFlat::shoot(engine, -M_PI, M_PI);
                    double xv = lxy * std::sin(rot);
                    double yv = lxy * std::cos(rot);
                    double zv = CLHEP::RandGaussQ::shoot(engine, 0., 5.);  // using gaussian for vz
                    double tv = CLHEP::RandGaussQ::shoot(engine, 0., 5.);  // using gaussian for ctau

                    vertex = new HepMC::GenVertex(HepMC::FourVector(xv, yv, zv, tv));
                    auto vertex_pos = vertex->point3d();

                    // Calculate eta
                    eta = calcEta(
                            charge / pt, eta_st2, phi, 
                            vertex_pos.x(), vertex_pos.y(), vertex_pos.z()
                    ); 
                }
            } else if (vertex_spectrum_ == "flatD0") { 
                // Calculate qinvpt
                double qinvpt = charge / pt;

                // Try to find vertex until eta_st2 is found
                double eta_st2 = eta;
                eta = 0;

                for (int itry = 0; (itry < 100) && (eta == 0); ++itry) {
                    // Delete Vertex
                    if (vertex != nullptr) {
                        delete vertex; 
                    }

                    // Calculate dxy
                    double rand_val = CLHEP::RandFlat::shoot(engine, 0., 1.);
                    double dxy = min_dxy_ + rand_val * (max_dxy_ - min_dxy_);

                    if (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5) {
                        dxy = -dxy;
                    }

                    // Calculate vertex
                    vertex = new HepMC::GenVertex(newVertexWithD0(engine, qinvpt, dxy, phi));
                    auto vertex_pos = vertex->point3d();

                    // Calculate eta
                    eta = calcEta(
                            qinvpt, eta_st2, phi, 
                            vertex_pos.x(), vertex_pos.y(), vertex_pos.z()
                    ); 
                }
            } else if (vertex_spectrum_ == "flatOneOverD0") { 
                // Calculate qinvpt
                double qinvpt = charge / pt;

                // Try to find vertex until eta_st2 is found
                double eta_st2 = eta;
                eta = 0;

                for (int itry = 0; (itry < 100) && (eta == 0); ++itry) {
                    // Delete Vertex
                    if (vertex != nullptr) {
                        delete vertex; 
                    }

                    // Generate dxy
                    double rand_val = CLHEP::RandFlat::shoot(engine, 0., 1.);
                    double dxy = min_one_over_dxy_ + rand_val * (max_one_over_dxy_ - min_one_over_dxy_);

                    if (CLHEP::RandFlat::shoot(engine, 0., 1.) < 0.5) {
                        dxy = -dxy;
                    }

                    if (dxy != 0.) {
                        dxy = 1. / dxy;
                    }

                    // Calculate vertex
                    vertex = new HepMC::GenVertex(newVertexWithD0(engine, qinvpt, dxy, phi));
                    auto vertex_pos = vertex->point3d();

                    // Calculate eta
                    eta = calcEta(
                            qinvpt, eta_st2, phi, 
                            vertex_pos.x(), vertex_pos.y(), vertex_pos.z()
                    ); 
                }
            } else {
                vertex = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0.));
            }
        }


        // Calculate data
        const HepPDT::ParticleData* particle_data = pdg_table.particle(
                HepPDT::ParticleID(abs(part_id_))
        );
        double mass = particle_data->mass().value();
        double theta = 2. * std::atan(std::exp(-eta));
        double mom = pt / std::sin(theta);
        double px = pt * std::cos(phi);
        double py = pt * std::sin(phi);
        double pz = pt / std::tan(theta);
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
        std::cout << " FlatRandomPtGunProducer2 : Event Generation Done " << std::endl;
    }
}

double FlatRandomPtGunProducer2::calcEta(
        const double& qinvpt,
        const double& eta_st2,
        const double& phi,
        const double& xv,
        const double& yv,
        const double& zv
) const {
    // Extract endcap
    double endcap = std::copysign(1, eta_st2);

    double eta_i = endcap * 2.0;
    double eta_im1 = 0;
    double y_i = 0;
    double y_im1 = 0;
    double init_step = 0.001 * endcap;
    double min_delta = 0.005;

    int    i_step = 0;
    bool   retry_evt = false;

    for (i_step = 0; i_step < 100; ++i_step) {
        // Calculate Eta at Station 2 for this step
        double eta_st2_i = calcEtaSt2(qinvpt, eta_i, phi, xv, yv, zv);

        // The problem to be minimized is the difference between eta_st2
        y_i = eta_st2 - eta_st2_i;

        // Short-Circuit: The error is below the threshold, consider this a good solution
        if (abs(y_i) < min_delta) {
            break;
        }

        // Short-Circuit: In order to perform optimization one step is required.
        // Use the first step to define a gradient for the next step.
        if (i_step == 0) {
            eta_im1 = eta_i;
            y_im1 = y_i;
            eta_i += init_step;
            continue;
        }

        // Calculate the gradient of y
        double dy = y_i - y_im1;
        double dx = eta_i - eta_im1;

        // Short-Circuit: There was no change in x, therefore this should be a minima
        if (dx == 0) {
            break;
        }

        // Current values are now previous ones
        eta_im1 = eta_i;
        y_im1 = y_i;

        // Calculate next eta value to test
        double dy_dx = dy / dx;
        eta_i = eta_i - y_i / dy_dx;

        // Short-Circuit: Solution beyon eta of 5, this is invalid 
        if (5 < abs(eta_i)) {
            retry_evt = true;
            break;
        }

        // Short-Circuit: Solution requires oposite signed eta, this is invalid
        if ((eta_i * endcap) < 0) {
            retry_evt = true;
            break;
        }
    }

    // Short-Circuit: Failed to find eta; return invalid eta
    if ((abs(y_i) > min_delta) and (i_step == 99)) {
        retry_evt = true;
    }

    // Short-Circuit: Eta is beyond solution region; return invalid eta
    if (retry_evt) {
        return 0;
    }

    return eta_i;
}

double FlatRandomPtGunProducer2::calcEtaSt2(
        const double& invpt,
        const double& eta,
        const double& phi,
        const double& x0,
        const double& y0,
        const double& z0
) const {
    // Constants
    const double B = 3.811;            // in Tesla, no unit conversion needed
    double ZST2 = 850. * cm;
    double Z4T = 650. * cm;

    if (eta < 0) {
        ZST2 *= -1;
        Z4T *= -1;
    } 

    // Propagate eta to station 2
    double rg = -1.0 / (0.003 * B * invpt) * cm; // rg = -pT/(0.003 q B)  [cm], radius of the circle
    double cot = std::sinh(eta);  // cot(theta), which is pz/pt
                                  
    double arg_term_4T, sin_term_4T, cos_term_4T;
    double arg_term_0T, sin_term_0T, cos_term_0T;

    if (std::abs(Z4T) < std::abs(ZST2)) {
        arg_term_4T = std::abs((Z4T - z0) / cot);  // with magfield
        sin_term_4T = (2 * rg) * std::sin(arg_term_4T / (2 * rg));  // with magfield
        cos_term_4T = (2 * rg) * (1 - std::cos(arg_term_4T / (2 * rg)));  // with magfield
        arg_term_0T = std::abs((ZST2 - Z4T) / cot);  // without magfield
        sin_term_0T = arg_term_0T;  // without magfield
        cos_term_0T = 0;  // without magfield
    } else {
        // Also need to check for the boundary at r where 4T -> 0T, ignore for now
        arg_term_4T = std::abs((ZST2 - z0) / cot);  // with magfield
        sin_term_4T = (2 * rg) * std::sin(arg_term_4T / (2 * rg));  // with magfield
        cos_term_4T = (2 * rg) * (1 - std::cos(arg_term_4T / (2 * rg)));  // with magfield
        arg_term_0T = 0;  // without magfield
        sin_term_0T = 0;  // without magfield
        cos_term_0T = 0;  // without magfield
    }

    double phi_st2 = phi + arg_term_4T / (2 * rg);  // phi at the boundary where 4T -> 0T
    double x_st2 = x0 + std::cos(phi) * sin_term_4T - std::sin(phi) * cos_term_4T + \
            std::cos(phi_st2) * sin_term_0T - std::sin(phi_st2) * cos_term_0T;
    double y_st2 = y0 + std::sin(phi) * sin_term_4T + std::cos(phi) * cos_term_4T + \
            std::sin(phi_st2) * sin_term_0T + std::cos(phi_st2) * cos_term_0T;
    double r_st2 = std::hypot(x_st2, y_st2);
    double cot_st2 = ZST2 / r_st2;
    double eta_st2 = std::asinh(cot_st2);

    return eta_st2;
}

HepMC::FourVector FlatRandomPtGunProducer2::newVertexWithD0(
        CLHEP::HepRandomEngine* engine,
        const double& invpt,
        const double& dxy,
        const double& phi
) const {
    // Constants
    const double B = 3.811;            // in Tesla, no unit conversion needed
    const double SIGMA_Z = 5.0 * cm; 

    // Declare variables
    double lxy, delta_ang; 
    double d0, alpha; 
    double d0_chk, lxy_chk;
    double xv, yv, zv, tv;

    // Extract dxy sign
    double dxy_sign = std::copysign(1, dxy);
    double udxy = std::abs(dxy);
        
    // Calculate the radius of gyration
    double rg = -1.0 / (0.003 * B * invpt) * cm;  // R = -pT/(0.003 q B)  [cm]
    double rg_sign = std::copysign(1, rg);
    double urg = std::abs(rg);

    // Calculate max/min lxy values
    double min_lxy, max_lxy;

    if (rg_sign == dxy_sign) {
        max_lxy = std::abs(udxy - 2 * urg);
    } else {
        max_lxy = std::abs(udxy + 2 * urg);
    }

    max_lxy = std::min(max_lxy, max_lxy_);
    min_lxy = udxy;
    
    // Try until lxy is less than max_lxy_
    int num_tries = 0;

    for (
            num_tries = 0; 
            num_tries == 0 || ((num_tries < 1000) && (max_lxy_ < lxy_chk)); 
            ++num_tries
    ) {
        // Generate flat lxy distribution
        lxy = CLHEP::RandFlat::shoot(engine, min_lxy, max_lxy);

        if (max_lxy_ < std::fabs(lxy)) {
            std::cout << "MuonGunLxyDebug lxy=" << (lxy / cm) 
                << " min lxy=" << (min_lxy / cm)
                << " max lxy=" << (max_lxy / cm)
                << " tries=" << num_tries
                << std::endl;
        }

        // Calculate delta angle
        double delta_ang_num = std::pow(lxy, 2) - std::pow(rg, 2) - std::pow(dxy - rg, 2); 
        double delta_ang_den = 2 * rg * (dxy - rg);

        delta_ang = std::acos(delta_ang_num / delta_ang_den);

        if (std::isnan(delta_ang)) {
            std::cout << "MuonGunAngDebug lxy=" << (lxy / cm) 
                << " min lxy=" << (min_lxy / cm)
                << " max lxy=" << (max_lxy / cm)
                << " delta num=" << (delta_ang_num / std::pow(cm, 2))
                << " delta den=" << (delta_ang_den / std::pow(cm, 2))
                << " tries=" << num_tries
                << std::endl;
        }

        // Calculate gyration center
        // Solve for the the cosine of the angle delta, 
        // since the cosine function is even there are two options.
        // 1. delta = alpha - phi -> alpha = delta + phi
        // 2. delta = phi - alpha -> alpha = phi - delta
        // For each attempt pick an option at random, they're both valid.
        double rand_val = CLHEP::RandFlat::shoot(engine, 0, 1);

        if (rand_val < 0.5) {
            alpha = phi - delta_ang; 
        } else {
            alpha = phi + delta_ang;
        }

        // Constrain alpha to [-pi, pi] range
        if (alpha < -M_PI) {
            alpha += (M_PI * 2);
        } else if (alpha >= M_PI) {
            alpha -= (M_PI * 2);
        }

        // Find the center of the gyration
        double xc = (dxy - rg) * std::sin(alpha);
        double yc = - (dxy - rg) * std::cos(alpha);

        // Find the production vertex with phi at the vertex
        xv = xc + rg * std::sin(phi);  // xc = xv - R sin(phi), note that xv is aX
        yv = yc - rg * std::cos(phi);  // yc = yv + R cos(phi), note that yv is aY
        
        // Check lxy
        lxy_chk = std::hypot(xv, yv);

        if (std::fabs(lxy - lxy_chk) >= 1) {
            std::cout << "MuonGunLxyMismatch lxy=" << (lxy / cm) 
                << " lxy_chk=" << (lxy_chk / cm)
                << " invpt=" << invpt
                << " rg=" << (rg / cm)
                << " d0=" << (d0 / cm)
                << " phi=" << phi
                << " alpha=" << alpha
                << " delta_ang=" << delta_ang
                << " min lxy=" << (min_lxy / cm)
                << " max lxy=" << (max_lxy / cm)
                << " tries=" << num_tries
                << std::endl;
        }

        // Calculate d0
        double xd0 = dxy * std::sin(alpha);
        double yd0 = - dxy * std::cos(alpha);

        double lxy_dot_d0 = xv * xd0 + yv * yd0;
        double lxy_d0_ang = std::acos(lxy_dot_d0 / lxy / udxy);

        if (lxy_d0_ang == 0) {
            d0 = 0;
        } else {
            double d0_sign;

            if (lxy_d0_ang < (M_PI / 2)) {
                d0_sign = 1;
            } else {
                d0_sign = -1;
            }

            d0 = d0_sign * udxy;
        }
        
        // Check d0
        double xc_test = xv - (rg * std::sin(phi));
        double yc_test = yv + (rg * std::cos(phi));

        double cg_test = std::hypot(xc_test, yc_test);
        double xc_norm = xc_test / cg_test;
        double yc_norm = yc_test / cg_test;

        double cg_minus_rg_test = cg_test - rg_sign * rg;
        double xd0_test = cg_minus_rg_test * xc_norm;
        double yd0_test = cg_minus_rg_test * yc_norm;
        double udxy_test = std::hypot(xd0_test, yd0_test);
    
        double lxy_dot_d0_test = xv * xd0_test + yv * yd0_test;
        double lxy_d0_ang_test = std::acos(lxy_dot_d0_test / lxy / udxy_test);

        if (lxy_d0_ang_test == 0) {
            d0_chk = 0;
        } else {
            double d0_chk_sign;

            if (lxy_d0_ang_test < (M_PI / 2)) {
                d0_chk_sign = 1;
            } else {
                d0_chk_sign = -1;
            }

            d0_chk = d0_chk_sign * udxy_test;
        }

        if (std::fabs(d0 - d0_chk) >= 1) {
            std::cout << "MuonGunD0Mismatch d0=" << (d0 / cm) 
                << " d0_chk=" << (d0_chk / cm)
                << " cg_test=" << (cg_test / cm)
                << " invpt=" << invpt
                << " rg=" << (rg / cm)
                << " lxy=" << (lxy / cm)
                << " phi=" << phi
                << " alpha=" << alpha
                << " delta_ang=" << delta_ang
                << " min lxy=" << (min_lxy / cm)
                << " max lxy=" << (max_lxy / cm)
                << " tries=" << num_tries
                << std::endl;
        }
    }
    
    if (num_tries > 1) {
        std::cout << "MuonGunAttempt invpt=" << invpt
            << " rg=" << (rg / cm)
            << " d0=" << (d0 / cm)
            << " lxy=" << (lxy / cm)
            << " phi=" << phi
            << " alpha=" << alpha
            << " delta_ang=" << delta_ang
            << " min lxy=" << (min_lxy / cm)
            << " max lxy=" << (max_lxy / cm)
            << " tries=" << num_tries
            << std::endl;
    }
    
    zv = CLHEP::RandGaussQ::shoot(engine, 0., SIGMA_Z);  // using gaussian for vz
    tv = CLHEP::RandGaussQ::shoot(engine, 0., SIGMA_Z);  // using gaussian for ctau

    return HepMC::FourVector(xv, yv, zv, tv);
}


//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(FlatRandomPtGunProducer2);
