#ifndef HAMILTONIAN_INL
#define HAMILTONIAN_INL

#include "model.inl"

template<class ParameterType>
class Hamiltonian
{
static_assert(std::is_base_of<ParameterBase, ParameterType>::value, "Parameter must inherit from ParameterSet for hamiltonians");
protected:
  const Model<ParameterType>& model;

  virtual void IntegratePath(
    ParameterType& parameters,
    ParameterType& momenta) const = 0;

public:
  explicit Hamiltonian(const Model<ParameterType>& model) : model(model) {}
  virtual ~Hamiltonian() {}

  void GenerateStep(ParameterType& parameter, ParameterType& momentum) const
  {
    IntegratePath(parameter, momentum);
  }

  virtual ParameterType RandomMomentum(const ParameterType& parameter, double temperature = 1.) const = 0;

  virtual double Energy(const ParameterType& parameter, const ParameterType& momentum) const = 0;
};

#endif
