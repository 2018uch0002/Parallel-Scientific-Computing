#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include "getCPU.h"

// generate random number from a uniform distribution
class RandomNumberGenrator{
    public:
        RandomNumberGenrator() {
            engine = std::mt19937(static_cast<unsigned int>(3564225786894));
            uni_distri = std::uniform_real_distribution<double>(0., 1.);
        }

        double rng() {
            return uni_distri(engine);
        }
    private:
        std::mt19937 engine; 
        std::uniform_real_distribution<double> uni_distri;

};

struct Particle{
    double x;
    double mu;
    bool is_alive;
};

class Tally{
    public:
        Tally (double xmin, double xmax, std::size_t N)
            : flux_gen_(), flux_avg_(), xmin_(xmin), xmax_(xmax), inv_dx_(), N_gen_(0) {
                flux_gen_ = std::vector<double>(N, 0.);
                flux_avg_ = std::vector<double>(N, 0.);
                inv_dx_ = static_cast<double>(N) / (xmax_ - xmin_);
            }

        void score_collision(const double x){
            std::size_t nx = static_cast<std::size_t>((x - xmin_) * inv_dx_);
            flux_gen_[nx] += 1. * inv_net_particles_;
        }

        void record_generation(){
            N_gen_++; 
            const double d_gen = static_cast<double>(N_gen_);
            for (std::size_t i = 0; i < flux_avg_.size(); i++){
                const double avg = flux_avg_[i];
                flux_avg_[i] = (avg * (d_gen-1.) + flux_gen_[i]) / d_gen;
            }
        }

        std::vector<double> flux() const { return flux_avg_; }
        std::size_t N_bins() const { return flux_avg_.size(); }
        double x_min() const { return xmin_; }
        double x_max() const { return xmax_; }

        void save_flux(std::string fname){
            std::ofstream file(fname);
            const double dx = 1./ inv_dx_;
            for (std::size_t i = 0; i < flux_avg_.size(); i++){
                file << xmin_ + (0.5 + static_cast<double>(i)) * dx << "\t" << flux_avg_[i] << "\n";
            }
            file.close();
        }

        void save_gen_flux(std::string fname){
            std::ofstream file(fname);
            const double dx = 1./ inv_dx_;
            for (std::size_t i = 0; i < flux_avg_.size(); i++){
                file << xmin_ + (0.5 + static_cast<double>(i)) * dx << "\t" << flux_gen_[i] << "\n";
            }
            file.close();
        }

        void gen_clear(){
            // fill zeros
            for (std::size_t i = 0; i < flux_avg_.size(); i++){
                flux_gen_[i] = 0.;
            }
        }

        void set_net_particles(std::size_t N){
            inv_net_particles_ = 1./ static_cast<double>(N);
        }

    private:
        std::vector<double> flux_gen_, flux_avg_;
        double xmin_, xmax_, inv_dx_, inv_net_particles_;
        std::size_t N_gen_;
};

int main(){
    double cpu0 = getCPU();
    // std::vector<double> Ea_values = {0.01, 0.03, 0.1};
    std::vector<double> Ea_values = {0.2, 0.3};
    std::vector<double> kavg_all = {0., 0., 0.};
    std::vector<double> kerr_all = {0., 0., 0.};
    for (std::size_t itr_Ea = 0; itr_Ea < 1; itr_Ea++){

        // random number generator
        RandomNumberGenrator rng;

        // cross section data
        const double Et = 1.;
        const double Ea = Ea_values[itr_Ea];
        const double Ef = 0.8 * Ea;
        const double nu = 2.5;
        const double Es = Et - Ea;

        const double probability_capture = Ea / Et;
        const double probability_scatter = 1. - Ea / Et;

        // total generations
        const std::size_t N_generations = 100;
        // total number of particles in one fisson generations
        std::size_t N_particles = 120000;

        // geometry
        const double x_min = 0.;
        const double x_max = 15.;

        // record the keff after the ceration generations
        const std::size_t inactive_gen = 20;
        double k_avg = 0.;
        double k_var = 0.;

        std::vector<double> k_eff(N_generations, 1.);
        Tally flux_estimate(x_min, x_max, 150);

        // store all particles
        std::vector<Particle> particle_bank;
        particle_bank.reserve(N_particles);
        for(std::size_t i = 0; i < N_particles; i++){
            Particle p;
            p.is_alive = true;
            p.x = rng.rng() * (x_max - x_min) + x_min;
            p.mu = -1. + 2. * rng.rng(); // born isotropically
            particle_bank.push_back(p);
        }

        
        // to store a new fission particle for net generation
        std::vector<Particle> next_gen_particle_bank; 


        // loop over each generations
        std::cout << "gen" << "\t" << "k" << "\t" << "k avg +/- var" << std::endl;
        for (std::size_t n_gen = 0; n_gen < N_generations; n_gen++){

            N_particles = particle_bank.size();
            double N_fission_particle = 0;
            double k_old = 1.;
            if (n_gen > 0) {
                k_old = k_eff[n_gen-1];
            }
            // flux_estimate.set_net_particles(N_particles);

            // loop over all particles in particular generation
            for (std::size_t i = 0; i < N_particles; i++){
                auto& p = particle_bank[i];

                while(p.is_alive == true){

                    // score collision
                    // flux_estimate.score_collision(p.x);
                    

                    const double s  = - log(rng.rng()) / Et;
                    p.x += s * p.mu;

                    // check if particle is in the domain
                    if (p.x < x_min || p.x > x_max ){
                        // particle leaked out
                        p.is_alive = false;
                        break;
                    } 

                    // check weather absorbed or leaked
                    const double xi = rng.rng();
                    
                    if (xi < probability_scatter){
                        // then scatter
                        // this problem is involves the isotropic
                        p.mu = -1. + 2. * rng.rng(); 
                    } else {
                        // absorbed
                        N_fission_particle += nu * Ef / Ea;
                        std::size_t I = std::floor(nu * Ef / (k_old * Ea) + rng.rng());
                        std::size_t t = 0;
                        while (t < I){
                            Particle p_new;
                            p_new.is_alive = true;
                            p_new.x = p.x;
                            p_new.mu = -1. + 2. * rng.rng();

                            next_gen_particle_bank.push_back(p_new);
                            t++;
                        }

                        p.is_alive = false;
                        break;
                    }
                }
            }

            // get the keff and store it
            double k_new = N_fission_particle / N_particles;
            k_eff[n_gen] = k_new;
            std::cout << n_gen + 1 << "\t" << k_new;

            // get the k_avg
            if (n_gen+1 > inactive_gen){
                const double delta_k = k_new - k_avg;
                const double d_gen = static_cast<double>(n_gen + 1 - inactive_gen);
                k_avg = ( k_new + k_avg * (d_gen - 1.) ) / d_gen;
                
                if (d_gen > 1.){
                    k_var = k_var + (delta_k * delta_k) / d_gen - k_var / (d_gen-1.);
                    std::cout << "\t" << k_avg << " +/- " << std::sqrt(k_var) << "\t N = " << N_particles <<"\tN_new =" << N_fission_particle;
                }

                // flux_estimate.record_generation();
                
            }
            // flux_estimate.gen_clear();

            std::cout << std::endl;
            particle_bank = next_gen_particle_bank;
            next_gen_particle_bank.clear();
            N_particles = particle_bank.size();
        }

        // double temp = 0.;
        // double temp_var = 0.;
        // for (std::size_t i = 20; i < N_generations; i++){
        //     temp += k_eff[i];
        // }
        // temp /= 80;
        // for (std::size_t i = 20; i < N_generations; i++){
        //     temp_var += (k_eff[i] - temp)*(k_eff[i] - temp);
        // }
        
        // std::cout << "\n\n" << temp<< "\t" << temp_var / 79 <<  std::endl;

        // flux_estimate.save_flux("flux_final_" + std::to_string(itr_Ea) + ".txt" ); 

        std::cout << "\nk_avg +/ err = " << k_avg << " +/- " << std::sqrt(k_var)  << std::endl;
        kavg_all[itr_Ea] = k_avg;
        kerr_all[itr_Ea] = std::sqrt(k_var);

        std::cout << "\n\n===================================\n\n" << std::endl;
    }

    std::cout << "TIME " << getCPU()-cpu0 << std::endl;
    std::ofstream file("keff.txt");
    for (std::size_t i = 0; i < Ea_values.size(); i++){
        file << Ea_values[i] << "\t" << kavg_all[i] << "\t" << kerr_all[i] << "\n";
    }
    file.close();

}