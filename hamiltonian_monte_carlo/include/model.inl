#ifndef MODEL_INL
#define MODEL_INL

#include "parameter_set.inl"
#include <iostream>

template<class ParameterType>
class Model
{
static_assert(std::is_base_of<ParameterBase, ParameterType>::value, "Parameter type must inherit ParameterSet for Model");
private:
  mutable int energy_evaluations = 0;
  mutable int partial_evaluations = 0;
  const double numerical_interval = 1e-7;

  virtual double CalculateEnergy(const ParameterType& parameters) const = 0;
  virtual ParameterType CalculateEnergyPartials(const ParameterType& parameters) const = 0;

  virtual ParameterType RealMapPartials(
      const ParameterType& mapped_parameters,
      const ParameterType& partials) const = 0;
public:
  explicit Model() {}
  virtual ~Model() {}

  typedef ParameterType parameter_type;
  virtual ParameterType ParameterMapReals(const ParameterType& parameter) const = 0;

  double Energy(const ParameterType& parameters) const 
  {
    energy_evaluations++;
    return CalculateEnergy(
        ParameterMapReals(parameters));
  }

  ParameterType EnergyPartials(const ParameterType& parameters) const
  {
    partial_evaluations++;
    ParameterType mapped_params = ParameterMapReals(
      parameters);

    return RealMapPartials(
        mapped_params,
        CalculateEnergyPartials(mapped_params));
  }

  int getEnergyEvaluationCount() const
  {
    return energy_evaluations;
  }

  int getEnergyPartialsEvaluationCount() const
  {
    return partial_evaluations;
  }

  virtual ParameterType getRandomInitialState() const = 0;

protected:

  ParameterType NumericalPartialEstimate(const ParameterType& parameter) const
  { 
    parameter_type perturbed_parameter;
    parameter_type partial_estimates;

    for (int i = 0; i < parameter_type::dimension; i++)
    {
      perturbed_parameter = parameter;

      perturbed_parameter.parameters.at(i) -= 0.5 * numerical_interval;
      const double old_energy =
        CalculateEnergy(perturbed_parameter);

      perturbed_parameter.parameters.at(i) += numerical_interval;
      const double new_energy =
        CalculateEnergy(perturbed_parameter);

      partial_estimates.parameters.at(i) =
        (new_energy - old_energy) / numerical_interval;

    }

    return partial_estimates;
  }
};

#endif
