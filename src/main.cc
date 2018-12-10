
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
#include <iomanip>
#include <fstream>

#include <random>
#include <string>

ExoplanetModel::parameter_type solar_model(
    std::string file_name,
    const ExoplanetModel::parameter_type& mass_matrix,
    double step_length = .1, 
    int path_length = 15, 
    int iterations = 2000, 
    int calibration_iterations = 20000, 
    int burn_in = 200,
    double burn_in_step = 20,
    int thinning_factor = 2,
    double tempering_factor = 1.,
    bool verbose = true,
    unsigned int repeats = 1)
{
  std::ofstream outfile("solar_model_sim.csv");
  std::ofstream calib_outfile("solar_model_sim_calib.csv");

  ExoplanetModel::parameter_type average_ess;
  std::vector<ExoplanetModel::parameter_type> means;
  std::vector<ExoplanetModel::parameter_type> variances;
  double average_acceptance_ratio = 0.;
  for (int rep = 0; rep < repeats; rep++)
  {
    std::cout << "Starting rep " << rep << std::endl;
    ExoplanetModel model(file_name);

    TrajectoryTemperedHamiltonian<ExoplanetModel::parameter_type>
      hamiltonian(model, step_length, path_length, tempering_factor);
    hamiltonian.setParameterMasses(mass_matrix);

    NUTS<ExoplanetModel::parameter_type>
      hmc(model, hamiltonian);

    SimpleHamiltonian<ExoplanetModel::parameter_type>
      burn_in_hamiltonian(model, burn_in_step, path_length);
    burn_in_hamiltonian.setParameterMasses(mass_matrix);

    NUTS<ExoplanetModel::parameter_type>
      burn_in_hmc(model, burn_in_hamiltonian);

    ExoplanetModel::parameter_type initial_state =
      model.getRandomInitialState();
    std::cout << "Initial state: " << initial_state << std::endl;

    std::cout << "Start time: " << model.getStartTime() << std::endl;
    std::cout << "Starting initial relaxation" << std::endl;
    //Move initial state to something fairly stable
    {
      burn_in_hmc.SimulateNSteps(burn_in, initial_state, verbose);
      burn_in_hmc.ClearHistory(); 
    }

    std::cout << "Starting calibration for ESS calculation" << std::endl;
    hmc.SimulateNSteps(calibration_iterations, initial_state, verbose);
    std::vector<ExoplanetModel::parameter_type> calibration =
      hmc.getSimulatedParameters(thinning_factor);
    hmc.ClearHistory();

    for (ExoplanetModel::parameter_type state : calibration)
      calib_outfile << std::setprecision(13) << model.ParameterMapReals(state) << std::endl;

    std::cout << "Running main simulation" << std::endl;
    hmc.SimulateNSteps(iterations, initial_state, verbose);
    //std::vector<ExoplanetModel::parameter_type> results =
    //  hmc.getSimulatedParameters(thinning_factor);

    for (ExoplanetModel::parameter_type state : hmc.getSimulatedParameters())
      outfile << std::setprecision(13) << model.ParameterMapReals(state) << std::endl;

    std::cout << "Calculating statistics" << std::endl;
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

    means.push_back(mean);
    variances.push_back(variance);

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

  std::cout << "Results: " << std::setprecision(6) << std::endl;
  if (verbose)
    std::cout << average_ess << "," << average_acceptance_ratio 
      << "," << repeats << "," << grand_mean << "," << grand_variance 
      << "," << pooled_variance << std::endl;

  return grand_mean;
}

template<int particles, class mcmc_type>
double lj_model(
    std::ofstream &outfile,
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
  for (int rep = 0; rep < repeats; rep++)
  {
    std::cout << "Starting rep " << rep << " with temp " << temperature << std::endl;
    LennardJonesModel<particles> model(x_size, y_size);

    TrajectoryTemperedHamiltonian<typename LennardJonesModel<particles>::parameter_type>
      hamiltonian(model, step_length, path_length, tempering_factor);

    mcmc_type hmc(model, hamiltonian, temperature);

    typename LennardJonesModel<particles>::parameter_type initial_state =
      model.getRandomInitialState();

    if(verbose)
      std::cout << "Relaxing state" << std::endl;
    //Move initial state to something fairly stable
    {
      SimpleHamiltonian<typename LennardJonesModel<particles>::parameter_type>
        hamiltonian(model, 0.0001, 100);
      HamiltonianMonteCarlo<typename LennardJonesModel<particles>::parameter_type>
        hmc(model, hamiltonian, temperature);

      hmc.SimulateNSteps(100, initial_state, verbose);
      hmc.ClearHistory(); 
    }

    if(verbose)
      std::cout << "Calibrating" << std::endl;
    hmc.SimulateNSteps(calibration_iterations, initial_state, verbose);
    std::vector<typename LennardJonesModel<particles>::parameter_type> calibration =
      hmc.getSimulatedParameters(thinning_factor);
    hmc.ClearHistory();

    std::vector<double> energy_calibration, energy_results;
    for (const typename LennardJonesModel<particles>::parameter_type& state : calibration)
    {
      energy_calibration.push_back(
        model.Energy(state));
    }

    if(verbose)
      std::cout << "Running main simulation" << std::endl;
    hmc.SimulateNSteps(iterations, initial_state, verbose);
    std::vector<typename LennardJonesModel<particles>::parameter_type> results =
      hmc.getSimulatedParameters(thinning_factor);

    for(const typename LennardJonesModel<particles>::parameter_type& state : results)
    {
      energy_results.push_back(
          model.Energy(state));
    }

    double effective_sample_size =
      MCMCDiagnostics::EffectiveSampleSize(energy_results, energy_calibration);

    double energy = MCMCDiagnostics::Mean(energy_results);

    std::vector<double> displacements =
      LennardJonesModel<particles>::CalculateMeanSquaredDisplacement(results);

    outfile << temperature << "," << energy 
    << "," << effective_sample_size << "," << hmc.getAcceptanceRatio()
    << "," << displacements.back() << std::endl;
  }
  return 0.;
}

int main()
{
  /*std::cout << "Starting Look Ahead" << std::endl;
  {
  std::ofstream outfile("lennard_jones_out_look_ahead_std_temp_long.csv");
  outfile << "T,E,ESS,Acceptance Ratio,Displacement Squared,Samples" << std::endl;
  for (int i = 1; i < 60; i++)
  {
  //step length, path length, iterations,
  //calibration iterations, thinning factor,
  //temperature, tempering factor, x size, y size,
  //verbose, repeats
    lj_model<16, LookAheadSampler<LennardJonesModel<16>::parameter_type> >(
        outfile,.001,10,16000,500,2,i * 0.25,1.,4.,4.,false, 30.);
  }
  }*/

  /*std::cout << "Starting Look Ahead with tempering" << std::endl;
  {
  std::ofstream outfile("lennard_jones_out_look_ahead_tempered_low_temp_long.csv");
  outfile << "T,E,ESS,Acceptance Ratio,Displacement Squared,Samples" << std::endl;
  for (int i = 1; i < 60; i++)
  {
  //step length, path length, iterations,
  //calibration iterations, thinning factor,
  //temperature, tempering factor, x size, y size,
  //verbose, repeats
    lj_model<16, LookAheadSampler<LennardJonesModel<16>::parameter_type> >(
        outfile,.001,10,16000,500,2,i * 0.05,1.05,4.,4.,false, 30.);
  }
  }*/

  std::cout << "Starting No U Turn" << std::endl;
  {
  std::ofstream outfile("lennard_jones_out_no_u_low_temp_long.csv");
  outfile << "T,E,ESS,Acceptance Ratio,Displacement Squared,Samples" << std::endl;
  for (int i = 1; i < 60; i++)
  {
  //step length, path length, iterations,
  //calibration iterations, thinning factor,
  //temperature, tempering factor, x size, y size,
  //verbose, repeats
    lj_model<16, NUTS<LennardJonesModel<16>::parameter_type> >(
        outfile,.025,1,16000,500,2,i * 0.05,1.,4.,4.,false, 30.);
  }
  }
  
  /*std::string file_name = "51_pegasi_256.txt";
  ExoplanetModel::parameter_type mass_matrix;
  mass_matrix.parameters.at(0) = 40;
  mass_matrix.parameters.at(1) = 1000;
  mass_matrix.parameters.at(2) = 20;
  mass_matrix.parameters.at(3) = 0.5;
  mass_matrix.parameters.at(4) = 10;
  mass_matrix.parameters.at(5) = 0.05;

  //file name, step length, path length, iterations
  //calibration iterations, burn in, burn in step, thinning factor, tempering factor
  //verbose, repeats
  solar_model(
      file_name, mass_matrix, .0000025, 8, 
      20000, 4000, 200, 0.0000025, 1,
      1.1, true, 1);*/

  
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
