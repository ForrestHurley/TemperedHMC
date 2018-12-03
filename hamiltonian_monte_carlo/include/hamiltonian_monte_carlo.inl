#ifndef HAMILTONIAN_MONTE_CARLO_INL
#define HAMILTONIAN_MONTE_CARLO_INL

#include "hamiltonian.inl"
#include "markov_chain_monte_carlo.inl"

template<class ParameterType>
class HamiltonianMonteCarlo : public MarkovChainMonteCarlo<ParameterType>
{
private:
  double temperature;
protected:
  const Hamiltonian<ParameterType>& hamiltonian;

  unsigned int reject_count = 0;
  unsigned int accept_count = 0;
public:
  HamiltonianMonteCarlo(
      const Model<ParameterType>& model, 
      const Hamiltonian<ParameterType>& hamiltonian, 
      double temperature = 1.) 
    : hamiltonian(hamiltonian),
      MarkovChainMonteCarlo<ParameterType>(model),
      temperature(temperature) {}

  virtual void SimulateStep(ParameterType& parameter) override
  {
    ParameterType old_parameter = parameter;
    ParameterType old_momentum = hamiltonian.RandomMomentum(parameter, temperature);
    ParameterType momentum = old_momentum;

    hamiltonian.GenerateStep(parameter, momentum);
    
    const double old_energy = hamiltonian.Energy(old_parameter, old_momentum);
    const double new_energy = hamiltonian.Energy(parameter, momentum);

    double probability =
      exp( (old_energy - new_energy ) / temperature); 

    if(this->getRandomUniform() < probability)
    {
      accept_count++;
      return;
    }

    reject_count++;
    parameter = old_parameter;
    return;
  }

  double getTemperature() const { return temperature; }
  void setTemperature(double new_temp) const { temperature = new_temp; }

  double getAcceptanceRatio() const 
  { return (double) accept_count / ( accept_count + reject_count ); }

  void ClearHistory()
  {
    MarkovChainMonteCarlo<ParameterType>::ClearHistory();
    reject_count = accept_count = 0;
  }
};

#endif
