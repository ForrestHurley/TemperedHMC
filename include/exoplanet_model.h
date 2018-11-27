#ifndef EXOPLANET_MODEL_H
#define EXOPLANET_MODEL_H

#include <vector>
#include <pair>
#include <string>

//Parameters
//Semi-major axis, eccentricity, argument of periapsis
//time of periastron passage, variance of measurements
class ExoplanetParameters : public ParameterSet<5>
{
public:
  double& getSemiMajorAxis() { return parameters.at(0); }
  void setSemiMajorAxis(double new_axis) { parameters.at(0) = new_axis; }

  double& getEccentricity() { return parameters.at(1); }
  void setEccentricity(double new_ecc) { parameters.at(1) = new_ecc; }

  double& getPeriapsis() { return parameters.at(2); }
  void setPeriapsis(double new_periapsis) { parameters.at(2) = new_periapsis; }
  
  double& getPeriastron() { return parameters.at(3); }
  void setPeriastron(double new_periastron) { parameters.at(3) = new_periastron; }

  double& getVariance() { return parameters.at(4); }
  void setVariance(double new_variance) { parameters.at(4) = new_variance; }
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
};

#endif 
