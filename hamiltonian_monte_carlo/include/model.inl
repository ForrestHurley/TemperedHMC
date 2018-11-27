#include "parameter_set.inl"

template<class ParameterType>
class Model<ParameterType>
{
static_assert(std::is_base_of<ParameterSet, ParameterType>::value, "Parameter type must inherit ParameterSet for Model");
private:
  int energy_evaluations = 0;
  int partial_evaluations = 0;

  virtual double CalculateEnergy(const ParameterType& parameters) const = 0;
  virtual ParameterType CalculateEnergyPartials(const ParameterType& parameters) const = 0;
public:
  explicit Model() {}
  virtual ~MarkovChainMonteCarlo() {}

  typedef ParameterType parameter_type;

  double Energy(const ParameterType& parameters) const 
  {
    energy_evaluations++;
    return CalculateEnergy(parameters);
  }

  ParameterType EnergyPartials(const ParameterType& parameters) const
  {
    partial_evaluations++;
    return CalculateEnergyPartials(parameters);
  }

  int getEnergyEvaluationCount() const
  {
    return energy_evaluations;
  }

  int getEnergyPartialsEvaluationCount() const
  {
    return partial_evalutations;
  }
};
