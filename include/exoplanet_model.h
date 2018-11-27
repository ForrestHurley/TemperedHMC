#ifndef EXOPLANET_MODEL_H
#define EXOPLANET_MODEL_H

#include <vector>
#include <pair>
#include <string>

//Parameters
//Semi-major axis, eccentricity, argument of periapsis
//time of periastron passage, variance of measurements
class ExoplanetParameters : public ParameterSet<6>
{
public:
  double& getSemiMajorAxis() { return parameters.at(0); }
  void setSemiMajorAxis(double new_axis) { parameters.at(0) = new_axis; }

  double& getEccentricity() { return parameters.at(1); }
  void setEccentricity(double new_ecc) { parameters.at(1) = new_ecc; }

  double& getPeriapsisLongitude() { return parameters.at(2); }
  void setPeriapsisLongitude(double new_periapsis) { parameters.at(2) = new_periapsis; }
  
  double& getPeriapsisTime() { return parameters.at(3); }
  void setPeriapsisTime(double new_periapsis) { parameters.at(3) = new_periapsis; }

  double& getVariance() { return parameters.at(4); }
  void setVariance(double new_variance) { parameters.at(4) = new_variance; }

  double& getPeriod() { return parameters.at(5); }
  void setPeriod(double new_period) { parameters.at(5) = new_period; }
};

class ExoplanetModel : public Model<ExoplanetParameters>
{
private:
  struct RadialVelocity
  {
    double velocity;
    double time;
  };

  std::vector<RadialVelocity> data_points;

public:
  explicit ExoplanetModel(std::vector<std::vector<double> > data_points);
  explicit ExoplanetModel(std::string filename);

private:
  double CalculateEnergy(const parameter_type& parameter) const override;
  parameter_type CalculateEnergyPartials(const parameter_type& parameter) const override;

  double ExpectedVelocity(const parameter_type& parameter, double time) const;
  parameter_type ExpectedVelocityPartials(const parameter_type& parameter, double time) const;

  double MassRatio(const parameter_type& parameter) const;
};

#endif 
