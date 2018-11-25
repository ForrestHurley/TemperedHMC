#ifndef MARKOV_CHAIN_MONTE_CARLO_INL
#define MARKOV_CHAIN_MONTE_CARLO_INL

#include "model.inl"

template<class ParamaterType>
class MarkovChainMonteCarlo<ParameterType>
{
static_assert(std::is_base_of<ParameterSet, ParameterType>::value, "Model parameter must inherit ParameterSet for MCMC");
private:
  std::vector<ParameterType> parameter_history;
protected:
  const Model& model;

public:
  explicit MarkovChainMonteCarlo(const Model& model) : model(model) {}
  virtual ~MarkovChainMonteCarlo() {}

  virtual void SimulateStep(ParameterType& parameter) = 0;

  void SimulateNSteps(unsigned int steps, ParameterType parameter)
  {
    parameter_history.clear();
    parameter_history.reserve(steps);

    for (int i = 0; i < steps; i++)
    {
      SimulateStep(parameter);
      parameter_history.push_back(parameter);
    }
  }

  const std::vector<ParameterType>&
    getSimulatedParameters() const { return parameter_history; }

  std::vector<ParameterType>
    getSimulatedParameters(int thinning_factor) const
  {
    std::vector<ParameterType> out;
    out.reserve(parameter_history.size() / thinning_factor);

    int remainder = parameter_history.size() % thinning_factor;
    
    for (int i = remainder; i < parameter_history.size(); i += thinning_factor)
      out.push_back(parameter_history.at(i));
    
    return out;
  }
};

#endif
