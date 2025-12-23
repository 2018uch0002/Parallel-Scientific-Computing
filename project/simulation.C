#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "simulation.h"
#include "getCPU.h"
#include "parseCommand.h"


// to define the threadBlockSize externally in .cu and update the value
// from *.C file
// ideal threadBlockSize should be either 128 or 256
int threadBlockSize = 256;

extern "C"
{  
  // function to generate the source particle for the first generation
  void generate_source_particles(Particle* &particle_bank, int& num_particles, 
          double& xmin, double& xmax, int& myrank, const int& n_bins, 
          double* &flux_hist);

  // Run k-eigenvalue simulation on GPU
  void run_gpu_simulation(Particle* &particles_bank, int& num_particles, double Et,
          double Ea, double Ef, double nu, double x_min, double x_max, 
          double& total_fission, int rank, int generation, double k_old, 
          const int inactive_gen, const int n_bins, double* &flux_hist);

  // clean the memory for particle-bank
  void free_particles(Particle* &particles);
  void free_flux(double* &flux);
}

// Helper functions
// Prints output for a flux histogram
void print_flux(double* &flux, int n_bins, double x_min, double dx) {
  for (int i = 0; i < n_bins; i++) {
    printf("Bin %d (< %lf): %lf\n", i, x_min + (0.5 + i) * dx, flux[i]);
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  string line;
  int total_N_particles = 1200;
  for(int n_arg=1; n_arg < argc; n_arg++){
      line = argv[n_arg];
      if(parseCommand(line, "-CUDAthread=", threadBlockSize)){}
      if(parseCommand(line, "-Nparticles=", total_N_particles)){}
  }

  // printf("PASS MPI_INIT\n"); //(DEBUG)

  // custom type because BYTES kept breaking
  // MPI_Datatype PARTICLE_MPI_TYPE;
  // {
  //   int blocklengths[3] = {1, 1, 1};
  //   MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  //   MPI_Aint displacements[3];

  //   displacements[0] = offsetof(Particle, x);
  //   displacements[1] = offsetof(Particle, mu);
  //   displacements[2] = offsetof(Particle, is_alive);

  //   MPI_Type_create_struct(3, blocklengths, displacements, types,
  //                          &PARTICLE_MPI_TYPE);
  //   MPI_Type_commit(&PARTICLE_MPI_TYPE);
  // }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Material Properties
  double Et = 1.0, Ea = 0.2, Ef = 0.16, nu = 2.5;

  // Geometry Details
  double x_min = 0.0, x_max = 15.0;

  // k-eigen simulation parameters
  const int total_gen = 100, inactive_gen = 20;
  const int active_gen = total_gen - inactive_gen;

  if (rank == 0)
    printf("Total Number of Paticles: %d\n", total_N_particles); //(DEBUG)

  // Flux Parameters
  const int n_bins = 20;
  const double dx = (x_max - x_min) / n_bins;

  // get the local_count for each mpi-rank
  int local_count = total_N_particles / size;
  if (rank == 0) {
    // add the extra particles at first rank
    local_count += total_N_particles % size;
  }

  if (rank == 0) {
    printf("====================================================\n");
    printf("k-eignvalue Monte Carlo Simulation on GPU\n\n");
    printf("   Generating the fission source particles.\n");
  }

  double *flux_avg = new double [n_bins]; 
  double *flux_hist;

  // generate the source particles uniformly and isotropically
  Particle *particles_bank;
  generate_source_particles(particles_bank, local_count, x_min, x_max, rank,
                            n_bins, flux_hist);

  MPI_Barrier(MPI_COMM_WORLD);

  // k-avg evaluation parameter
  double *k_generation;
  if (rank == 0) {
    k_generation = new double[active_gen];
  }
  double k_new = 1.;
  double k_avg = 0.0, k_var = 0.0;

  uint64_t start_cycle, total_cycle = 0;

  if (rank == 0) {
    printf("Particles: %d, Inactive Gen: %d, Active Gen: %d\n\n",
           total_N_particles, inactive_gen, active_gen);
    printf("Gen.\t k-gen\t kavg +/- std\n");
  }

  double cpu0 = getCPU();
  for (int gen = 0; gen < total_gen; gen++) {
    start_cycle = getCPU();

    // get the total particles
    int global_particles;
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Allreduce(&local_count, &global_particles, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    

    double local_fission = 0;
    run_gpu_simulation(particles_bank, local_count, Et, Ea, Ef, nu, x_min,
                       x_max, local_fission, rank, gen, k_new, inactive_gen,
                       n_bins, flux_hist);

    MPI_Barrier(MPI_COMM_WORLD);

    // Get the total average fission particle
    double global_fission;
    double *flux_global;
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Allreduce(&local_fission, &global_fission, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    
    if ((gen > inactive_gen) && (size > 1)) {
      MPI_Reduce(flux_hist, flux_global, n_bins, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
    } // Combine the flux calculation histograms from all ranks
    MPI_Barrier(MPI_COMM_WORLD);
    
    k_new = global_fission / ((double)global_particles);

    if (rank == 0) {
      // // << FLUX-CAL  Calculate kth average flux
      // if (gen > inactive_gen) {
      //   for (int i = 0; i < n_bins; i++) {
      //     flux_avg[i] =
      //         (flux_avg[i] * (gen - inactive_gen - 1) + flux_global[i]) /
      //         (gen - inactive_gen);
      //   }
      // }

      // Output stats
      printf("%d\t %2.5lf\t", gen, k_new);
      k_generation[gen] = k_new;

      if (gen > inactive_gen) {
        double delta_k = k_new - k_avg;
        const double d_gen = (double)(gen - inactive_gen);
        k_avg = (k_new + k_avg * (d_gen - 1.)) / d_gen;

        if (d_gen > 1.) {
          k_var = k_var + (delta_k * delta_k) / d_gen - k_var / (d_gen - 1.);
          printf("%2.5lf +/- %2.5lf", k_avg, sqrt(k_var));
        }
      }
      printf("\n");m
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }

  double total_time = getCPU() - cpu0;
  // final output
  if (rank == 0) {
    printf("\n\n===============================================\n");
    printf("Final k_avg: %.5f +/- %.5f\n", k_avg, sqrt(k_var));
    printf("Total simulation time: %.2f seconds\n", total_time);
    printf("===============================================\n\n");

    // print_flux(flux_avg, n_bins, x_min, dx); // << FLUX-CAL
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // free(particles); // we cannot free it here we have to free on cuda
  free_particles(particles_bank);
  free_flux(flux_hist);
  // free(flux_avg);

  MPI_Finalize();

  return 0;
}
