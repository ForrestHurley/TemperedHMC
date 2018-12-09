
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "exoplanet_model.h"

TEST(ExoplanetDerivativeTest, Zero)
{
  std::vector<std::vector<double> > test_points;
  test_points.push_back(std::vector<double>{0, 0});

  ExoplanetModel test_model(test_points);

  ExoplanetModel::parameter_type test_parameter;

  ExoplanetModel::parameter_type energy_partials =
    test_model.EnergyPartials(test_parameter);
  
  std::cout << energy_partials << std::endl;

  ExoplanetModel::parameter_type partial_estimates;

  double interval = 1e-7;
  for (int i = 0; i < test_parameter.dimension; i++)
  {
    ExoplanetModel::parameter_type perturbed_parameter =
      test_parameter;

    perturbed_parameter.parameters.at(i) -= 0.5 * interval;
    const double old_energy =
      test_model.Energy(perturbed_parameter);

    perturbed_parameter.parameters.at(i) += interval;
    const double new_energy =
      test_model.Energy(perturbed_parameter);

    partial_estimates.parameters.at(i) =
      (new_energy - old_energy) / interval;

    EXPECT_NEAR(
        partial_estimates.parameters.at(i),
        energy_partials.parameters.at(i),
        1e-7);
  }

}

TEST(ExoplanetDerivativeTest, HalfValueParam)
{
  std::vector<std::vector<double> > test_points;
  test_points.push_back(std::vector<double>{0, 0});

  ExoplanetModel test_model(test_points);

  ExoplanetModel::parameter_type test_parameter(2.);

  ExoplanetModel::parameter_type energy_partials =
    test_model.EnergyPartials(test_parameter);

  std::cout << energy_partials << std::endl;

  ExoplanetModel::parameter_type partial_estimates;

  double interval = 1e-7;
  for (int i = 0; i < test_parameter.dimension; i++)
  {
    ExoplanetModel::parameter_type perturbed_parameter =
      test_parameter;

    perturbed_parameter.parameters.at(i) -= 0.5 * interval;
    const double old_energy =
      test_model.Energy(perturbed_parameter);

    perturbed_parameter.parameters.at(i) += interval;
    const double new_energy =
      test_model.Energy(perturbed_parameter);

    partial_estimates.parameters.at(i) =
      (new_energy - old_energy) / interval;

    EXPECT_NEAR(
        partial_estimates.parameters.at(i),
        energy_partials.parameters.at(i),
        1e-5);
  }

}

TEST(ExoplanetDerivativeTest, NinetyDegree)
{

}
