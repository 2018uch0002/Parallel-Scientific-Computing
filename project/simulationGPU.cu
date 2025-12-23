#include "simulation.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/scan.h>


// extern theredBlockSize
extern int threadBlockSize;

extern "C"{
  
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

// function to generate the source to begin the generations for k-eignvalue simulation
__global__ 
void generate_source_particles_kernel(Particle* particle_bank, 
        int num_particles, double xmin, double xmax, int myrank){
  
  // Generate the particles which are born isotropically and uniformly
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  curandState local_state;
  __uint64_t seed = ((__uint64_t)(1<<31) - 1) - (((__uint64_t)(1<<20) - 1)) - (1023*myrank);
  curand_init(seed, idx, 0, &local_state);

  if (idx < num_particles){
    Particle p;
    // uniform distributed between xmin and xmax
    p.x = xmin + (xmax - xmin) * curand_uniform_double(&local_state);
    // born isotropically
    p.mu = -1. + 2. * curand_uniform_double(&local_state);
    p.is_alive = 1;
    particle_bank[idx] = p;
  }
}

void generate_source_particles(Particle* &particle_bank, int& num_particles, 
          double& xmin, double& xmax, int& myrank, const int& n_bins, 
          double* &flux_hist){
  // get the cuda device
  int device_count;
  cudaGetDeviceCount(&device_count);

  // set the cuda device
  cudaSetDevice(myrank % device_count);

  // Allocate memory for the particle bank
  cudaMallocManaged((void **)& (particle_bank), num_particles * sizeof(Particle));
  
  // Allocate memory for the flux histogram
  cudaMallocManaged((void**)& (flux_hist), n_bins * sizeof(double));
  cudaDeviceSynchronize();

  int num_block = (num_particles + threadBlockSize -1) / threadBlockSize ;

  generate_source_particles_kernel<<<num_block, threadBlockSize>>>(particle_bank, 
                                                                    num_particles, 
                                                                    xmin, xmax, myrank);
  cudaDeviceSynchronize();
}


__global__ 
void neutron_transport_kernel(Particle *particles, int num_particles,
                                         double Et, double Ea, double Ef,
                                         double nu, double x_min, double x_max,
                                         double probability_scatter, int *new_counts,
                                         double *d_total_fission,
                                         double k_old, int gen, int myrank, const int inactive_gen, 
                                         const int n_bins, double* flux_hist) {
  
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_particles)
    return;

  curandState local_state;
  __uint64_t seed = ((__uint64_t)(1<<31) - 1) - (((__uint64_t)(1<<21) - 1) * 2 * (gen + 1)) - (1023*myrank);
  curand_init(seed, idx, 0, &local_state);

  // curandState local_state = states[idx];
  Particle& p = particles[idx];
  new_counts[idx] = 0;

  while (p.is_alive) {
    // score the collision estimator
    //Update Flux Calculation Histogram
    if(gen > inactive_gen){
      atomicAdd(flux_hist + int((p.x - x_min) * (n_bins / (x_max - x_min))), 1);
    } 

    // particle transport
    double s = -log(curand_uniform_double(&local_state)) / Et;
    p.x += s * p.mu;

    // check the boundary condition
    // only vacuum BC is applied
    if (p.x < x_min || p.x > x_max) {
      p.is_alive = 0;
      break;
    }
  
    if (curand_uniform_double(&local_state) < probability_scatter) {
      // scatter the particle isotropically
      p.mu = -1.0 + 2.0 * curand_uniform_double(&local_state);
    } else {
      // absorption
      double num_avg_born = nu * Ef / (k_old * Ea);
      double avg_born_particle = nu * Ef / Ea;
      atomicAdd(d_total_fission, avg_born_particle);
      
      int produced = (int)(num_avg_born + curand_uniform_double(&local_state));
      new_counts[idx] = produced;
      p.is_alive = 0;
      break;

    }
  }
  // particles[idx] = p;
}

__global__ 
void generateNewParticles(Particle* new_bank, int *prefix_sum,
                                     Particle *old_particles, int *new_counts,
                                     int num_particles,
                                     double x_min, double x_max, int gen, int myrank) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= num_particles)
    return;

  int count = new_counts[idx];
  if (count == 0)
    return;

  // curandState local_state = states[idx];
  curandState local_state;
  __uint64_t seed = ((__uint64_t)(1<<31) - 1) - (((__uint64_t)(1<<21) - 1) * gen) - (1023*myrank);
  curand_init(seed, idx, 0, &local_state);
  int start_idx = prefix_sum[idx];

  for (int i = 0; i < count; i++) {
    Particle p;
    p.x = old_particles[idx].x;
    p.mu = -1.0 + 2.0 * curand_uniform_double(&local_state);
    p.is_alive = 1;
    new_bank[start_idx + i] = p;
  }
}

void run_gpu_simulation(Particle* &particles_bank, int& num_particles, double Et,
        double Ea, double Ef, double nu, double x_min, double x_max, 
        double& total_fission, int rank, int generation, double k_old, 
        const int inactive_gen, const int n_bins, double* &flux_hist) {

  // get the cuda device
  int device_count;
  cudaGetDeviceCount(&device_count);

  // set the cuda device
  cudaSetDevice(rank % device_count);

  // CUDA memory management
  int *d_prefix, *d_counts;
  double *d_fission;

  const int N = num_particles;
  const int threads = threadBlockSize;
  const int blocks = (N + threads - 1) / threads;

  // allocations on GPU
  cudaMalloc(&d_counts, N * sizeof(int));
  cudaMalloc(&d_prefix, (N + 1) * sizeof(int));
  // cudaMalloc(&d_states, N * sizeof(curandState));
  cudaMalloc(&d_fission, sizeof(double));


  // copy particles' count to device
  cudaMemset(d_counts, 0, N * sizeof(int));
  cudaMemset(d_fission, 0, sizeof(double));
  cudaMemset(flux_hist, 0, n_bins * sizeof(double));
  
  // run transport kernel
  double p_scatter = 1. - Ea/Et;
  neutron_transport_kernel<<<blocks, threads>>>(particles_bank, N, Et, Ea, Ef, nu,
                                                x_min, x_max, p_scatter,
                                                d_counts, d_fission, k_old, generation, 
                                                rank, inactive_gen, n_bins, flux_hist);
  cudaDeviceSynchronize();

  thrust::exclusive_scan(thrust::cuda::par.on(0), d_counts, d_counts + (N + 1), d_prefix);
  cudaDeviceSynchronize();

  // Particle for new fission generaitons
  int total_new;
  cudaMemcpy(&total_new, d_prefix + N, sizeof(int), cudaMemcpyDeviceToHost);

  Particle *new_particle_bank;
  cudaMallocManaged((void **)& new_particle_bank, total_new * sizeof(Particle));

  // new particles
  if (total_new > 0) {
    //cudaMalloc(&d_new_bank, total_new * sizeof(Particle));
    generateNewParticles<<<blocks, threads>>>(new_particle_bank, d_prefix, particles_bank, 
                                            d_counts, N, x_min, x_max, generation, rank);
  }

  cudaDeviceSynchronize();

  // copy results back
  cudaMemcpy(&total_fission, d_fission, sizeof(double),
             cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();

  // update particles
  cudaFree(particles_bank);
  num_particles = total_new;
  if (total_new > 0) {
    particles_bank = new_particle_bank;
    // cudaMemcpy(*particles, d_new_bank, total_new * sizeof(Particle),
    //            cudaMemcpyDeviceToHost);
  }

  cudaDeviceSynchronize();

  // FREE Memory
  cudaFree(d_counts);
  cudaFree(d_prefix);
  // cudaFree(d_states);
  cudaFree(d_fission);
}


void free_particles(Particle* &particles){
  cudaFree(particles);
}


void free_flux(double* &flux){
  cudaFree(flux);
}