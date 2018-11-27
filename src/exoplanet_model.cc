#include "exoplanet_model.h"

ExoplanetModel::ExoplanetModel(std::vector<std::vector<double> > data_points)
{

}

ExoplanetModel::ExoplanetModel(std::string filename)
{

}

double ExoplanetModel::CalculateEnergy(const ExoplanetModel::parameter_type& parameter) const
{
  double energy = 0.;
  
  for(RadialVelocity datum : data_points)
  {

  }

  return energy; 
}

ExoplanetModel::parameter_type ExoplanetModel::CalculateEnergyPartials(
  const ExoplanetModel::parameter_type& parameter) const
{
  parameter_type partials;

  for(radialVelocity datum : data_points)
  {

  }

  return partials;
}
