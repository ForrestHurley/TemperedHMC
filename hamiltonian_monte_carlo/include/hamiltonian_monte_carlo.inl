#ifndef HAMILTONIAN_MONTE_CARLO_INL
#define HAMILTONIAN_MONTE_CARLO_INL

#include "hamiltonian.inl"

template<class ParameterType>
class HamiltonianMonteCarlo<ParameterType> : public MarkovChainMonteCarlo<ParameterType>
{
protected:
  const Hamiltonian& hamiltonian;

  virtual double KineticEnergy(const ParameterType& momentum) const = 0;

public:
  HamiltonianMonteCarlo(const Model& model, const Hamiltonian& hamiltonian) 
    : hamiltonian(hamiltonian),
      MarkovChainMonteCarlo(model) {}

  virtual void SimulateStep(ParameterType& parameter) override
  {
    ParameterType old_parameter = parameter;
    ParameterType old_momentum = hamiltonian.RandomMomentum();
    ParameterType momentum = old_momentum;

    hamiltonian.GenerateStep(parameter, momentum)
    
    double probability =
      exp(model.Energy(old_parameter) - model.Energy(parameter) + 
          KineticEnergy(old_momentum) - KineticEnergy(momentum));

    if(getRandomUniform() < probability)
      return;

    parameter = old_parameter;
    return;
  }
};

#endif
