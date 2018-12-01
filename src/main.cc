
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

int main()
{
  //mean, variance
  GaussianModel model = GaussianModel(0., 1.);

  //model, step length, path length
  SimpleHamiltonian<GaussianModel::parameter_type> 
    hamiltonian(model, 1., 1);

  //model, hamiltonian, temperature
  HamiltonianMonteCarlo<GaussianModel::parameter_type>
    mcmc(model, hamiltonian, 1.);

  GaussianModel::parameter_type initial_state;

  mcmc.SimulateNSteps(10000, initial_state);
  std::vector<double> calibration = mcmc.getSimulatedParameters(1, 0);
  mcmc.ClearHistory();

  //steps, intial_state
  mcmc.SimulateNSteps(1000, initial_state);

  //thinning factor, parameter index
  std::vector<double> results = mcmc.getSimulatedParameters(1, 0);

  double mean = MCMCDiagnostics::Mean(results);
  double variance = MCMCDiagnostics::Variance(results, mean);

  double effective_sample_size =
    MCMCDiagnostics::EffectiveSampleSize(results, calibration);
  
  std::cout << "Mean: " << mean << " Variance: " << variance << std::endl;
  std::cout << "Effective Sample Size: " << effective_sample_size << std::endl;
  std::cout << "Acceptance Ratio: " << mcmc.getAcceptanceRatio() << std::endl; 

  return 0;
}
