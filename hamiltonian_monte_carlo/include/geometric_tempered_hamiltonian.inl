#ifndef GEOMETRIC_TEMPERED_HAMILTONIAN
#define GEOMETRIC_TEMPERED_HAMILTONIAN

#include "simple_hamiltonian.inl"

template<class ParameterType>
class GeometricTemperedHamiltonian<ParameterType> : public SimpleHamiltonian<ParameterType>
{
private:
  mutable ParameterType current_parameter; 
  //horrible hack, requires significant refactoring to fix cleanly (or code repetition)

  double temperature = 1;

public:
  GeometricTemperedHamiltonian(const Model& model, const Hamiltonian& hamiltonian)
    : SimpleHamiltonian(model, hamiltonian) {}

  virtual void IntegratePath(
    ParameterType parameter&,
    ParameterType momentum&) const override;

  virtual double Energy(const ParameterType& parameter, const ParameterType& momentum) const override;

  virtual ParameterType getParameterMasses() override { return G(current_parameter); }
  virtual ParameterType G(const ParameterType& parameter) const;

  virtual ParameterType RandomMomentum(const ParameterType& parameter) const override;

  void setTemperature(double new_temperature) { temperature = new_temperature; }
  double getTemperature() const { return temperature; }
};

template<class ParameterType>
void GeometricTemperedHamiltonian<ParameterType>::IntegratePath(ParameterType& parameter, ParameterType& momentum)
{

}

template<class ParameterType>
double GeometricTemperedHamiltonian<ParameterType>::Energy(const ParameterType& parameter, const ParameterType& momentum)
{
  current_parameter = parameter;
  const ParameterType dimension_energy = momentum * momentum / ( 2 * getParameterMasses());
  const double kinetic = dimension_energy.Sum();
  const double potential = model.Energy(parameter);
  const double mass_energy = log(sqrt(getParameterMasses().Magnitude()));
  return kinetic + mass_energy + potential;
}

template<class ParameterType>
ParameterType GeometricTemperedHamiltonian<ParameterType>::G(const ParameterType& parameter) const
{
  const double probability = exp(-model.Energy(parameter));
  const double mass = pow(probability, 1. - 1. / temperature);
  return ParameterType(mass);
}

template <class ParameterType>
ParameterType GeometricTemperedHamiltonain<ParameterType>::RandomMomentum(const ParameterType& parameter) const
{
  current_parameter = parameter;
  return SimpleHamiltonian<ParameterType>::RandomMomentum(parameter);
}


#endif
