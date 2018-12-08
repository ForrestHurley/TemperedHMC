
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
#include <string>

ExoplanetModel::parameter_type solar_model(
    std::string file_name,
    double step_length = .1, 
    int path_length = 15, 
    int iterations = 2000, 
    int calibration_iterations = 20000, 
    int thinning_factor = 2,
    double tempering_factor = 1.,
    bool verbose = true,
    unsigned int repeats = 1)
{
  ExoplanetModel::parameter_type average_ess;
  std::vector<ExoplanetModel::parameter_type> means;
  std::vector<ExoplanetModel::parameter_type> variances;
  double average_acceptance_ratio = 0.;
  for (int rep = 0; rep < repeats; rep++)
  {
    ExoplanetModel model(file_name);

    SimpleHamiltonian<ExoplanetModel::parameter_type>
      hamiltonian(model, step_length, path_length);

    NUTS<ExoplanetModel::parameter_type>
      hmc(model, hamiltonian);

    ExoplanetModel::parameter_type initial_state =
      model.getRandomInitialState();

    //Move initial state to something fairly stable
    {
      SimpleHamiltonian<ExoplanetModel::parameter_type>
        hamiltonian(model, 0.0001, 100);
      HamiltonianMonteCarlo<ExoplanetModel::parameter_type>
        hmc(model, hamiltonian);

      hmc.SimulateNSteps(100, initial_state);
      hmc.ClearHistory(); 
    }

    hmc.SimulateNSteps(calibration_iterations, initial_state);
    std::vector<ExoplanetModel::parameter_type> calibration =
      hmc.getSimulatedParameters(thinning_factor);
    hmc.ClearHistory();

    hmc.SimulateNSteps(iterations, initial_state);
    //std::vector<ExoplanetModel::parameter_type> results =
    //  hmc.getSimulatedParameters(thinning_factor);

    //loop over all parameters and summarize
    ExoplanetModel::parameter_type effective_sample_size;
    ExoplanetModel::parameter_type mean;
    ExoplanetModel::parameter_type variance;
    for (int i = 0; i < effective_sample_size.dimension; i++)
    {
      mean.parameters.at(i) =
        MCMCDiagnostics::Mean(
            hmc.getSimulatedParameters(thinning_factor, i));
      variance.parameters.at(i) =
        MCMCDiagnostics::Variance(
            hmc.getSimulatedParameters(thinning_factor, i),
            mean.parameters.at(i));
      effective_sample_size.parameters.at(i) = 
        MCMCDiagnostics::EffectiveSampleSize(
            hmc.getSimulatedParameters(thinning_factor, i),
            hmc.getSimulatedParameters(thinning_factor, i, calibration));
    }
    //double effective_sample_size =
    //  MCMCDiagnostics::EffectiveSampleSize();

    average_ess += effective_sample_size;
    average_acceptance_ratio += hmc.getAcceptanceRatio();

  }
  average_ess /= repeats;
  average_acceptance_ratio /= repeats;

  ExoplanetModel::parameter_type pooled_variance;
  for (const ExoplanetModel::parameter_type& state : variances)
    pooled_variance += state;
  pooled_variance /= variances.size();

  ExoplanetModel::parameter_type grand_mean;
  for (const ExoplanetModel::parameter_type& state : means)
    grand_mean += state;
  grand_mean /= means.size();

  ExoplanetModel::parameter_type grand_variance;
  for (const ExoplanetModel::parameter_type& state : means)
    grand_variance += (grand_mean - state) * (grand_mean - state);
  grand_variance /= means.size() - 1;

  if (verbose)
    std::cout << "," << average_ess << "," << average_acceptance_ratio 
      << "," << repeats << "," << grand_mean << "," << grand_variance 
      << "," << pooled_variance << std::endl;
  return grand_mean;
}

template<int particles>
double lj_model(
    double step_length = .1, 
    int path_length = 15, 
    int iterations = 2000, 
    int calibration_iterations = 20000, 
    int thinning_factor = 2,
    double temperature = 1.,
    double tempering_factor = 1.,
    double x_size = 16.,
    double y_size = 16.,
    bool verbose = true,
    unsigned int repeats = 1)
{
  double average_displacement = 0.;
  double average_energy = 0.;
  double average_ess = 0.;
  double average_acceptance_ratio = 0.;
  for (int rep = 0; rep < repeats; rep++)
  {
    LennardJonesModel<particles> model(x_size, y_size);

    SimpleHamiltonian<typename LennardJonesModel<particles>::parameter_type>
      hamiltonian(model, step_length, path_length);

    NUTS<typename LennardJonesModel<particles>::parameter_type>
      hmc(model, hamiltonian, temperature);

    typename LennardJonesModel<particles>::parameter_type initial_state =
      model.getRandomInitialState();

    //Move initial state to something fairly stable
    {
      SimpleHamiltonian<typename LennardJonesModel<particles>::parameter_type>
        hamiltonian(model, 0.0001, 100);
      HamiltonianMonteCarlo<typename LennardJonesModel<particles>::parameter_type>
        hmc(model, hamiltonian, temperature);

      hmc.SimulateNSteps(100, initial_state);
      hmc.ClearHistory(); 
    }

    hmc.SimulateNSteps(calibration_iterations, initial_state);
    std::vector<typename LennardJonesModel<particles>::parameter_type> calibration =
      hmc.getSimulatedParameters(thinning_factor);
    hmc.ClearHistory();

    std::vector<double> energy_calibration, energy_results;
    for (const typename LennardJonesModel<particles>::parameter_type& state : calibration)
    {
      energy_calibration.push_back(
        model.Energy(state));
    }

    hmc.SimulateNSteps(iterations, initial_state);
    std::vector<typename LennardJonesModel<particles>::parameter_type> results =
      hmc.getSimulatedParameters(thinning_factor);

    for(const typename LennardJonesModel<particles>::parameter_type& state : results)
    {
      energy_results.push_back(
          model.Energy(state));
    }

    double effective_sample_size =
      MCMCDiagnostics::EffectiveSampleSize(energy_results, energy_calibration);

    average_energy += MCMCDiagnostics::Mean(energy_results);

    std::vector<double> displacements =
      LennardJonesModel<particles>::CalculateMeanSquaredDisplacement(results);

    average_displacement += displacements.back();
    average_ess += effective_sample_size;
    average_acceptance_ratio += hmc.getAcceptanceRatio();

  }
  average_energy /= repeats;
  average_displacement /= repeats;
  average_ess /= repeats;
  average_acceptance_ratio /= repeats;

  if (verbose)
    std::cout << temperature << "," << average_energy 
    << "," << average_ess << "," << average_acceptance_ratio
    << "," << average_displacement << "," << repeats << std::endl;
  return average_energy;
}

int main()
{
  //std::cout << "T,E,ESS,Acceptance Ratio,Displacement Squared,Samples" << std::endl;
  /*for (int i = 1; i < 40; i++)
  {
    lj_model<16>(.025,1,500,500,2,i,1.,4.,4.,true, 1.);
  }*/

  
  std::string file_name = "exodata.txt";
  //file name, step length, path length, iterations
  //calibration iterations, thinning factor, tempering factor
  //verbose, repeats
  ExoplanetModel::parameter_type solar_model(
      file_name, .005, 30, 
      2000, 20000, 2,
      1., true, 1);

  
  /*//typedef LennardJonesModel<2> ModelType;
  typedef GaussianModel ModelType;
  //mean, variance
  GaussianModel model = GaussianModel(0., 1.);

  //LennardJonesModel<2> model;

  const double step_length = .1;
  const int path_length = 15;
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
  */

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
