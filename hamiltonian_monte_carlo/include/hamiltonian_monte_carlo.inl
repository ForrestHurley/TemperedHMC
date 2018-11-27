#ifndef HAMILTONIAN_MONTE_CARLO_INL
#define HAMILTONIAN_MONTE_CARLO_INL

#include "hamiltonian.inl"

template<class ParameterType>
class HamiltonianMonteCarlo<ParameterType> : public MarkovChainMonteCarlo<ParameterType>
{
protected:
  const Hamiltonian& hamiltonian;

public:
  HamiltonianMonteCarlo(const Model& model, const Hamiltonian& hamiltonian) 
    : hamiltonian(hamiltonian),
      MarkovChainMonteCarlo(model) {}

  virtual void SimulateStep(ParameterType& parameter) override
  {
    ParameterType old_parameter = parameter;
    ParameterType old_momentum = hamiltonian.RandomMomentum(parameter);
    ParameterType momentum = old_momentum;

    hamiltonian.GenerateStep(parameter, momentum)
    
    double probability =
      exp(hamiltonian.Energy(old_parameter, old_momentum) - 
          hamiltonian.Energy(parameter, momentum)); 

    if(getRandomUniform() < probability)
      return;

    parameter = old_parameter;
    return;
  }
};

#endif
