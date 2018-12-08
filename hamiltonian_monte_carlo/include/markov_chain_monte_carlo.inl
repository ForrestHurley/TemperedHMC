#ifndef MARKOV_CHAIN_MONTE_CARLO_INL
#define MARKOV_CHAIN_MONTE_CARLO_INL

#include "model.inl"
#include "parameter_set.inl"
#include <random>

template<class ParameterType>
class MarkovChainMonteCarlo
{
static_assert(std::is_base_of<ParameterBase, ParameterType>::value, "Model parameter must inherit ParameterSet for MCMC");
private:
  std::vector<ParameterType> parameter_history;
protected:
  const Model<ParameterType>& model;

  static double getRandomUniform()
  {
    static thread_local std::random_device device;
    static thread_local std::mt19937_64 generator(device());
    static thread_local std::uniform_real_distribution<double> uniform(0., 1.);

    return uniform(generator);
  }

public:
  explicit MarkovChainMonteCarlo(const Model<ParameterType>& model) : model(model) {}
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

  std::vector<double>
    getSimulatedParameters(int thinning_factor, int parameter_index) const
  {
    std::vector<double> out;
    out.reserve(parameter_history.size());

    int remainder = parameter_history.size() % thinning_factor;

    for (int i = remainder; i < parameter_history.size(); i += thinning_factor)
      out.push_back(parameter_history.at(i).parameters.at(parameter_index));
    
    return out;
  }

  static std::vector<double>
    getSimulatedParameters(
        int thinning_factor, int parameter_index,
        const std::vector<ParameterType>& external_data)
  {
    std::vector<double> out;
    out.reserve(external_data.size());

    int remainder = external_data.size() % thinning_factor;

    for (int i = remainder; i < external_data.size(); i+= thinning_factor)
      out.push_back(external_data.at(i).parameters.at(parameter_index));

    return out;
  }

  void ClearHistory()
  {
    parameter_history = std::vector<ParameterType>();
  }

};

#endif
