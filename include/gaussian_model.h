#ifndef GAUSSIAN_MODEL_H
#define GAUSSIAN_MODEL_H

#include "model.inl"

class GaussianModel : public Model<ParameterSet<1> >
{
private:
  double mean, variance;

public:
  GaussianModel(double mean, double variance) : 
    mean(mean), variance(variance) {}

  parameter_type ParameterMapReals(const parameter_type& parameter) const
  { return parameter; }

  parameter_type getRandomInitialState() const;

private:
  double CalculateEnergy(const parameter_type& parameters) const override;
  parameter_type CalculateEnergyPartials(
      const parameter_type& parameters) const override;

  parameter_type RealMapPartials(
      const parameter_type& mapped_parameters,
      const parameter_type& partials) const override
  { return partials; }
};

#endif

