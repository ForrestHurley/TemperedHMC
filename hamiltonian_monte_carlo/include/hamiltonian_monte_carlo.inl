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

  void SimulateStep(ParameterType& parameter) override
  {
    hamiltonian.GenerateStep(parameter)
  }
};

#endif
