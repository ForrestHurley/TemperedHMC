#include "exoplanet_model.h"
#include <cassert>
#include <iostream>
#include <fstream>

ExoplanetModel::ExoplanetModel(std::vector<std::vector<double> > data_points)
{
  for (std::vector<double> point : data_points)
  {
    assert(point.size() == 2, "Please ensure that the exoplanet model is being passed n x 2 data");
    this->data_points.push_back(
        RadialVelocity(point.at(0), point.at(1)));
  }
}

ExoplanetModel::ExoplanetModel(std::string filename)
{
  std::ifstream file(filename);
  std::string line;
  std::vector<std::vector<double> > data_points;

  while(std::getline(file, line))
  {
    std::vector<double> line_data;
    std::stringstream line_stream(line);

    double value;
    while(line_stream >> value)
    {
      line_data.push_back(value
    }
    if (line_data.size() != 2)
    {
      std::cerr << "Lines must have exactly two data points on them" << std:: endl;
      assert(false);
    }

    data_points.push_back(line_data);
  }

  for (std::vector<double> point : data_points)
  {
    this->data_points.push_back(
        RadialVelocity(point.at(0), point.at(1)));
  }
}

double ExoplanetModel::PlanetaryMass(
    const ExoplanetModel::parameter_type& real_parameter) const
{
  const parameter_type parameter = ParameterMapReals(real_parameter)
  const double stellar_mass = parameter.getStellarMass();
  const double inverse_n = parameter.getPeriod() / (2 * pi);
  const double gravitational_parameter = stellar_mass * G;
  const double power_ratio = 
    pow( gravitational_parameter / parameter.getSemiMajorAxis(), 3. / 2.);

  const double planet_mass = 
    inverse_n * power_ratio - stellar_mass;
  return planet_mass;
}

double ExoplanetModel::CalculateEnergy(const ExoplanetModel::parameter_type& parameter) const
{
  double energy = PriorEnergy(parameter);
  
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
  parameter_type partials = PriorEnergyPartials(parameter);

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
  const double x_speed =
    StellarXVelocity(parameter, time);
  const double y_speed =
    StellarYVelocity(parameter, time);

  const double radial_velocity =
    x_speed * cos(parameter.getPeriapsisLongitude()) +
    y_speed * sin(parameter.getPeriapsisLongitude());
  return radial_velocity;
}

ExoplanetModel::parameter_type ExoplanetModel::ExpectedVelocityPartials(
  const ExoplanetModel::parameter_type parameter, double time) const
{

  const double x_speed =
    StellarXVelocity(parameter, time);
  const double y_speed =
    StellarYVelocity(parameter, time);
  const double theta =
    parameter.getPeriapsisLongitude();

  parameter_type velocity_partials =
    cos(theta) * StellarXVelocityPartials(parameter, time) +
    sin(theta) * StellarYVelocityPartials(parameter, time);

  velocity_partials.getPeriapsisLongitude() =
    - x_speed * sin(theta) + y_speed * cos(theta);

  return velocity_partials;
}

double ExoplanetModel::EccentricAnomaly(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  const double n = 2 * pi / parameter.getPeriod();
  const double M = n * ( time - parameter.getPeriapsisTime() );
  const double eccentricity = parameter.GetEccentricity();

  double E = M;
  double E_old = E + 1;

  while (abs(E - E_old) > 1e-6)
  {
    E_old = E;
    
    //https://arxiv.org/pdf/1609.00915.pdf
    E = ( M + eccentricity * ( sin(E_old) - E_old * cos(E_old) ) ) /
      ( 1 - eccentricity * cos(E_old) );
  }

  return E;
}

ExoplanetModel::parameter_type ExoplanetModel::EccentricAnomalyPartials(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  parameter_type eccentric_partials;

  const double anomaly = EccentricAnomaly(parameter, time);

  eccentric_partials.getEccentricity() =
    sin(anomaly) / ( 1 - parameter.getEccentricity() * cos(anomaly) );

  eccentric_partials.getPeriod() =
    2 * pi * ( time - parameter.getPeriapsisTime() ) /
    ( parameter.getPeriod() * parameter.getPeriod() *
      ( parameter.getEccentricity() * cos(anomaly) - 1. ) );

  eccentric_partials.getPeriapsisTime() =
    2 * pi / ( parameter.getPeriod() *
        ( parameter.getEccentricity() * cos(anomaly) - 1. ) );

  return eccentric_partials;
}

double ExoplanetModel::StellarSemiMajorAxis(
    const ExoplanetModel::parameter_type& parameter) const
{
  const double two_body_sm_axis =
    G * parameter.getPeriod() * parameter.getStellarMass() / 
    ( 4 * pi * pi * parameter.getSemiMajorAxis() );

  //a1 + a2 = a
  const double stellar_sm = 
    sqrt(two_body_sm_axis) - parameter.getSemiMajorAxis();
  return stellar_sm;
}

ExoplanetModel::parameter_type ExoplanetModel::StellarSemiMajorAxisPartials(
    const ExoplanetModel::parameter_type& parameter) const
{
  parameter_type axis_partials;
  const double a = parameter.getSemiMajorAxis();
  const double p = parameter.getPeriod();

  const double gravitational_parameter = G * parameter.getStellarMass();
  const double ratio = gravitational_parameter / ( 16 * pi * pi * a );

  axis_partials.getPeriod() = sqrt(ratio / p ) - a;
  axis_partials.getSemiMajorAxis() = - sqrt(ratio * p / (a * a) ) - 1;

  return axis_partials;
}

double ExoplanetModel::StellarXVelocity(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  const double n = 2 * pi / parameter.getPeriod();
  const double anomaly = EccentricAnomaly(parameter, time);
  const double stellar_sm_axis = StellarSemiMajorAxis(parameter);

  const double x_speed =
    - stellar_sm_axis * n * sin(anomaly) / 
    ( 1 - parameter.getEccentricity() * cos(anomaly) );
  return x_speed;
}

ExoplanetModel::parameter_type ExoplanetModel::StellarXVelocityPartials(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  const double anomaly = EccentricAnomaly(parameter, time);
  const double stellar_axis = StellarSemiMajorAxis(parameter);
  const double sin_E = sin(anomaly);
  const double cos_E = cos(anomaly);
  const double e = parameter.getEccentricity();
  const double P = parameter.getPeriod();
  const double T = parameter.getPeriapsisTime();
  const double n = 2 * pi / P;
  const double common_denominator = 1. - e * cos_E;

  const parameter_type anomaly_partials = 
    EccentricAnomalyPartials(parameter, time);
  const parameter_type stellar_partials =
    StellarSemiMajorAxisPartials(parameter); 

  parameter_type x_velocity_partials;
  
  const double stellar_d_axis = stellar_partials.getSemiMajorAxis();
  x_velocity_partials.getSemiMajorAxis() =
    - n * stellar_d_axis * sin_E / common_denominator;

  const double stellar_d_period = stellar_partials.getPeriod();
  const double anomaly_d_period = anomaly_partials.getPeriod();
  x_velocity_partials.getPeriod() = 
    n * ( P * stellar_d_period * sin_E * ( - common_denominator ) + 
        stellar_axis * ( P * anomaly_d_period * (e - cos_E) + 
          sin_E * common_denominator)) /
    ( P * common_denominator * common_denominator );

  const double anomaly_d_ecc = anomaly_partials.getEccentricity();
  x_velocity_partials.getEccentricity() =
    n * stellar_axis * ( anomaly_d_ecc * ( e - cos_E ) - sin_E * cos_E ) /
    ( common_denominator * common_denominator );

  const double anomaly_d_peri_time = anomaly_partials.getPeriapsisTime();
  x_velocity_partials.getPeriapsisTime() =
    n * stellar_axis * anomaly_d_peri_time * ( e - cos_E ) /
    ( common_denominator * common_denominator );

  return x_velocity_partials;
}

double ExoplanetModel::StellarYVelocity(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  const double n = 2 * pi / parameter.getPeriod();
  const double anomaly = EccentricAnomaly(parameter, time);
  const double stellar_sm_axis = StellarSemiMajorAxis(parameter);
  const double e = parameter.getEccentricity();

  const double y_speed =
    - stellar_sm_axis * n * cos(anomaly) * sqrt( 1 - e * e ) / 
    ( 1 - parameter.getEccentricity() * cos(anomaly) );
  return y_speed;
}

ExoplanetModel::parameter_type ExoplanetModel::StellarYVelocityPartials(
    const ExoplanetModel::parameter_type& parameter, double time) const
{
  const double anomaly = EccentricAnomaly(parameter, time);
  const double stellar_axis = StellarSemiMajorAxis(parameter);
  const double sin_E = sin(anomaly);
  const double cos_E = cos(anomaly);
  const double e = parameter.getEccentricity();
  const double P = parameter.getPeriod();
  const double T = parameter.getPeriapsisTime();
  const double n = 2 * pi / P;
  const double common_denominator = 1. - e * cos_E;
  const double axes_ratio = sqrt( 1 - e * e );

  const parameter_type anomaly_partials = 
    EccentricAnomalyPartials(parameter, time);
  const parameter_type stellar_partials =
    StellarSemiMajorAxisPartials(parameter); 

  parameter_type y_velocity_partials;

  const double stellar_d_axis = stellar_partials.getSemiMajorAxis();
  y_velocity_partials.getSemiMajorAxis() =
    n * axes_ratio * stellar_d_axis * cos_E * P / common_denominator;

  const double stellar_d_period = stellar_partials.getPeriod();
  const double anomaly_d_period = anomaly_partials.getPeriod();
  y_velocity_partials.getPeriod() =
    - n * axes_ratio * ( P * stellar_d_period * cos_E * ( - common_denominator ) + 
        stellar_axis * ( P * anomaly_d_period * sin_E - cos_E * ( e * cos_E + 1 ) ) ) /
    ( P * common_denominator * common_denominator );

  const double anomaly_d_ecc = anomaly_partials.getEccentricity();
  y_velocity_partials.getEccentricity() =
    n * stellar_axis * ( ( e * e - 1 ) * anomaly_d_ecc * sin_E + 
        cos_E * ( cos_E - e ) ) /
    ( axes_ratio * common_denominator * common_denominator );

  const double anomaly_d_peri_time = anomaly_partials.getPeriapsisTime();
  y_velocity_parialts.getPeriapsisTime() = 
    - n * axes_ratio * stellar_axis * anomaly_d_peri_time * sin_E /
    ( common_denominator * common_denominator );

  return y_velocity_partials;
}


ExoplanetModel::parameter_type ExoplanetModel::ParameterMapReals(
    const parameter_type& parameter) const
{
  parameter_type out;

  out.setSemiMajorAxis
}

ExoplanetModel::parameter_type ExoplanetModel::RealMapParameters(
    const parameter_type& parameter) const
{

}


