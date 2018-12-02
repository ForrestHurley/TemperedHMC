
#include "model.inl"
#include "hamiltonian_monte_carlo.inl"
#include "no_u_turn_sampler.inl"
#include "look_ahead_sampler.inl"
#include "simple_hamiltonian.inl"
#include "trajectory_tempered_hamiltonian.inl"
#include "mcmc_diagnostics.h"

#include "gaussian_model.h"
#include "exoplanet_model.h"
#include "lennard_jones_model.inl"

#include <vector>
#include <iostream>

#include <random>

int main()
{
  typedef LennardJonesModel<2> ModelType;
  //mean, variance
  //GaussianModel model = GaussianModel(0., 1.);

  LennardJonesModel<2> model;

  const double step_length = 1.;
  const int path_length = 4;
  const int iterations = 2000;

  //model, step length, path steps
  SimpleHamiltonian<ModelType::parameter_type> 
    hamiltonian(model, step_length, path_length);
  SimpleHamiltonian<ModelType::parameter_type>
    calib_hamiltonian(model, 0.01, 100);

  //model, step length, path steps, alpha ratio
  TrajectoryTemperedHamiltonian<ModelType::parameter_type>
    tempered_hamiltitonian(model, step_length, path_length, 1.05);

  //model, hamiltonian, temperature
  HamiltonianMonteCarlo<ModelType::parameter_type>
    calib_mcmc(model, calib_hamiltonian, 1.);
  //HamiltonianMonteCarlo<GaussianModel::parameter_type>
  //  mcmc(model, hamiltonian, 1.);

  SimpleHamiltonian<ModelType::parameter_type>
    nuts_hamiltonian(model, step_length, 1);
  NUTS<ModelType::parameter_type>
    mcmc(model, nuts_hamiltonian);

  //model, hamiltonian, maximum steps
  //LookAheadSampler<GaussianModel::parameter_type>
  //  mcmc(model, hamiltonian, 5);

  ModelType::parameter_type initial_state
    = model.getRandomInitialState();

  std::cout << "Starting calibration iterations" << std::endl;

  calib_mcmc.SimulateNSteps(20000, initial_state);
  std::vector<double> calibration = calib_mcmc.getSimulatedParameters(2, 0);
  calib_mcmc.ClearHistory();

  //steps, intial_state
  mcmc.SimulateNSteps(iterations, initial_state);

  //thinning factor, parameter index
  std::vector<double> results = mcmc.getSimulatedParameters(2, 0);

  double mean = MCMCDiagnostics::Mean(results);
  double variance = MCMCDiagnostics::Variance(results, mean);

  double effective_sample_size =
    MCMCDiagnostics::EffectiveSampleSize(results, calibration);
  
  std::cout << "Mean: " << mean << " Variance: " << variance << std::endl;
  std::cout << "Effective Sample Size: " << effective_sample_size << std::endl;
  std::cout << "Total Iterations: " << iterations << std::endl;
  std::cout << "Acceptance Ratio: " << mcmc.getAcceptanceRatio() << std::endl; 
  

  /*std::random_device device;
  std::mt19937_64 twister(device());
  std::uniform_real_distribution<double> distribution(0., 4.);

  std::vector<double> calibration, results;
  for (int i = 0; i < 1000; i++)
  {
    calibration.push_back(i * 0.001 + 1.5);
    results.push_back(distribution(twister));
  }

  double effective_sample_size =
    MCMCDiagnostics::EffectiveSampleSize(results, calibration);

  std::cout << "ESS: " << effective_sample_size << std::endl;
  */
  return 0;
}
