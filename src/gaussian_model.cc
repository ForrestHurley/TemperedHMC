#include "gaussian_model.h"

#include <random>

GaussianModel::parameter_type GaussianModel::getRandomInitialState() const
{
  static thread_local std::random_device device;
  static thread_local std::mt19937_64 twister(device());
  static thread_local std::normal_distribution<double> normal(0., 1.);

  parameter_type out;
  out.parameters.at(0) = normal(twister);
  return out;
}

double GaussianModel::CalculateEnergy(
    const GaussianModel::parameter_type& parameters) const
{
  const double difference = parameters.parameters.at(0) - mean;
  return difference * difference / ( 2. * variance);
}

GaussianModel::parameter_type GaussianModel::CalculateEnergyPartials(
    const GaussianModel::parameter_type& parameters) const
{
  const double difference = parameters.parameters.at(0) - mean;
  parameter_type out;
  out.parameters.at(0) = difference / variance;
  return out;
}

