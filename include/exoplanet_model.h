#ifndef EXOPLANET_MODEL_H
#define EXOPLANET_MODEL_H

#include <vector>
#include <string>
#include "parameter_set.inl"
#include "model.inl"

#include "math.h"

//Parameters
//Semi-major axis, eccentricity, argument of periapsis
//time of periastron passage, variance of measurements
class ExoplanetParameters : public ParameterSet<6>
{
public:
  ExoplanetParameters() {}
  ExoplanetParameters(double val) : ParameterSet<6>(val) {}
  ExoplanetParameters(const ParameterSet<6>& old) : ParameterSet<6>(old) {}

  double& getSemiMajorAxis() { return parameters.at(0); }
  double getSemiMajorAxis() const { return parameters.at(0); }
  void setSemiMajorAxis(double new_axis) { parameters.at(0) = new_axis; }

  double& getEccentricity() { return parameters.at(1); }
  double getEccentricity() const { return parameters.at(1); }
  void setEccentricity(double new_ecc) { parameters.at(1) = new_ecc; }

  double& getPeriapsisLongitude() { return parameters.at(2); }
  double getPeriapsisLongitude() const { return parameters.at(2); }
  void setPeriapsisLongitude(double new_periapsis) { parameters.at(2) = new_periapsis; }

  double& getPeriapsisTime() { return parameters.at(3); }
  double getPeriapsisTime() const { return parameters.at(3); }
  void setPeriapsisTime(double new_periapsis) { parameters.at(3) = new_periapsis; }

  double& getVariance() { return parameters.at(4); }
  double getVariance() const { return parameters.at(4); }
  void setVariance(double new_variance) { parameters.at(4) = new_variance; }

  double& getPeriod() { return parameters.at(5); }
  double getPeriod() const { return parameters.at(5); }
  void setPeriod(double new_period) { parameters.at(5) = new_period; }

  double stellar_mass = 1;
  double& getStellarMass() { return stellar_mass; }
  double getStellarMass() const { return stellar_mass; }
  void setStellarMass(double new_mass) { stellar_mass = new_mass; }
};

class ExoplanetModel : public Model<ExoplanetParameters>
{
private:
  constexpr static double pi = 3.14159265358979323846264;
  constexpr static double G = 39.478;

  double start_time;

  struct RadialVelocity
  {
    double time;
    double velocity;

    RadialVelocity(double time, double velocity) :
      velocity(velocity), time(time) {}
  };

  std::vector<RadialVelocity> data_points;

public:
  explicit ExoplanetModel(std::vector<std::vector<double> > data_points);
  explicit ExoplanetModel(std::string filename);

  double PlanetaryMass(const parameter_type& parameter) const;
  virtual parameter_type ParameterMapReals(const parameter_type& parameter) const override;

  parameter_type getRandomInitialState() const override;

  double getStartTime() const { return start_time; }

private:
  double CalculateEnergy(const parameter_type& parameter) const override;
  parameter_type CalculateEnergyPartials(const parameter_type& parameter) const override;

  double ExpectedVelocity(const parameter_type& parameter, double time) const;
  parameter_type ExpectedVelocityPartials(const parameter_type& parameter, double time) const;

  double EccentricAnomaly(const parameter_type& parameter, double time) const;
  parameter_type EccentricAnomalyPartials(const parameter_type& parameter, double time) const;

  double StellarSemiMajorAxis(const parameter_type& parameter) const;
  parameter_type StellarSemiMajorAxisPartials(const parameter_type& parameter) const;

  double StellarXVelocity(const parameter_type& parameter, double time) const;
  parameter_type StellarXVelocityPartials(const parameter_type& parameter, double time) const;

  double StellarYVelocity(const parameter_type& parameter, double time) const;
  parameter_type StellarYVelocityPartials(const parameter_type& parameter, double time) const;

  double PriorEnergy(const parameter_type& parameter) const;
  parameter_type PriorEnergyPartials(const parameter_type& parameter) const;

  static double Logit(double in) { return log( in / ( 1. - in ) ); }
  static double InvLogit(double in) { return exp( in ) / ( 1. + exp( in ) ); }
  static double LogitDeriv(double in) { return 1. / ( in * ( 1. - in ) ); }

  virtual parameter_type RealMapPartials(
      const parameter_type& mapped_parameters,
      const parameter_type& partials) const override;

  //Friend classes for unit testing
  //friend class ExoplanetDerivativeTest_Zero_Test;
  //friend class ExoplanetDerivativeTest_NinetyDegree_Test;
};

#endif 
