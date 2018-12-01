#include "gaussian_model.h"

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

