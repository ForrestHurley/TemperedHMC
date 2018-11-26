#ifndef GEOMETRIC_TEMPERED_HAMILTONIAN
#define GEOMETRIC_TEMPERED_HAMILTONIAN

#include "simple_hamiltonian.inl"

template<class ParameterType>
class GeometricTemperedHamiltonian<ParameterType> : public SimpleHamiltonian<ParameterType>
{
public:
  GeometricTemperedHamiltonian(const Model& model, const Hamiltonian& hamiltonian)
    : SimpleHamiltonian(model, hamiltonian) {}

  virtual void SimulateStep(ParameterType parameter) override;
};

template<class ParameterType>
void GeometricTemperedHamiltonian<ParameterType>::SimulateStep(ParameterType parameter)
{

}


#endif
