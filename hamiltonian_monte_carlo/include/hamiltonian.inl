#ifndef HAMILTONIAN_INL
#define HAMILTONIAN_INL

#include "model.inl"

template<class ParameterType>
class Hamiltonian<ParameterType>
{
static_assert(std::is_base_of<ParameterSet, ParameterType>::value, "Parameter must inherit from ParameterSet for hamiltonians");
protected:
  const Model& model;

  virtual void IntegratePath(
    ParameterType parameters&,
    ParameterType momenta&) const = 0;
  virtual ParameterType RandomMomentum() const = 0;

public:
  explicit Hamiltonian(const Model& model) : model(model) {}
  virtual ~Hamiltonian() {}

  void GenerateStep(ParameterType& parameter) const
  {
    ParameterType momentum = RandomMomentum();
    IntegratePath(parameter, momentum);
  }

};

#endif
