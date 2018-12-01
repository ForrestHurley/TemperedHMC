#ifndef MODEL_INL
#define MODEL_INL

#include "parameter_set.inl"

template<class ParameterType>
class Model
{
static_assert(std::is_base_of<ParameterBase, ParameterType>::value, "Parameter type must inherit ParameterSet for Model");
private:
  mutable int energy_evaluations = 0;
  mutable int partial_evaluations = 0;

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
        parameters,
        CalculateEnergyPartials(parameters));
  }

  int getEnergyEvaluationCount() const
  {
    return energy_evaluations;
  }

  int getEnergyPartialsEvaluationCount() const
  {
    return partial_evaluations;
  }
};

#endif
