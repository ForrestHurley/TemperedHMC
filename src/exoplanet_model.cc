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
    const double expected_radial_velocity = 
      ExpectedVelocity(parameter, datum.time);
    
    //probability is proportional to
    //e ^ ( - (x - mu) ^2 / ( 2 sigma ^ 2 ))
    //energy is -log(probabilty) or
    // (x - mu) ^ 2 / (2 sigma ^ 2 )
    const double difference = expected_radial_velocity - datum.velocity;
    const double datum_energy =
      difference * difference / ( 2 * parameter.getVariance() );
      
    energy += datum_energy;
  }

  return energy; 
}

ExoplanetModel::parameter_type ExoplanetModel::CalculateEnergyPartials(
  const ExoplanetModel::parameter_type& parameter) const
{
  parameter_type partials;

  for(radialVelocity datum : data_points)
  {
    const parameter_type expected_velocity_partials =
      ExpectedVelocityPartials(parameter, datum.time);
    const double expected_radial_velocity =
      ExpectedVelocity(parameter, datum.time);
    const double difference = expected_radial_velocity - datum.velocity;
    
    //For everything except variance
    // d/d theta E = 2 (f(theta) - a) * d / d theta f / (2 variance)
    parameter_type datum_partials =
      difference * expected_velocity_partials / parameter.getVariance();

    datum_partials.getVariance() = 
      - difference * difference / (2 * parameter.getVariance() * parameter.getVariance());

    partials += datum_partials;
  }

  return partials;
}

double ExoplanetModel::ExpectedVelocity(
  const ExoplanetModel::parameter_type parameter, double time) const
{

}

ExoplanetModel::parameter_type ExoplanetModel::ExpectedVelocityPartials(
  const ExoplanetModel::parameter_type parameter, double time) const
{

}

double ExoplanetModel::MassRatio(const ExoplanetModel::parameter_type& parameter) const
{

}
