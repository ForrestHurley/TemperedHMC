#ifndef SIMPLE_HAMILTONIAN_INL
#define SIMPLE_HAMILTONIAN_INL

#include "hamiltonian.inl"
#include <random>

template<class ParameterType>
class SimpleHamiltonian : public Hamiltonian<ParameterType>
{
private:
  double step_length;
  unsigned int path_length;
  ParameterType parameter_masses;

public:
  explicit SimpleHamiltonian(
    const Model<ParameterType>& model,
    double step_length = 0.01,
    int path_length = 100)
      : Hamiltonian<ParameterType>(model),
        step_length(step_length),
        path_length(path_length), 
        parameter_masses(1.) { }

  virtual ~SimpleHamiltonian() {}  

  void setStepLength(double new_length) { step_length = new_length; }
  double getStepLength() const { return step_length; }

  void setPathLength(unsigned int new_length) { path_length = new_length; }
  unsigned int getPathLength() const { return path_length; }

  void setParameterMasses(const ParameterType& new_masses) 
    { parameter_masses = new_masses; }
  virtual const ParameterType& getParameterMasses() const { return parameter_masses; }

  virtual void IntegratePath(
    ParameterType& parameters,
    ParameterType& momenta) const override;

  void LeapfrogStep(
    ParameterType& parameters, 
    ParameterType& momenta) const;

  virtual ParameterType RandomMomentum(const ParameterType& parameter, double temperature = 1.) const override;

  virtual double Energy(
      const ParameterType& parameter, 
      const ParameterType& momentum) const override;
};

template<class ParameterType>
double SimpleHamiltonian<ParameterType>::Energy(
    const ParameterType& parameter, const ParameterType& momentum) const
{
  const ParameterType dimension_energy = momentum * momentum / ( getParameterMasses() * 2. );
  const double kinetic = dimension_energy.Sum();
  const double potential = this->model.Energy(parameter);
  return kinetic + potential;
}

template<class ParameterType>
void SimpleHamiltonian<ParameterType>::IntegratePath(
  ParameterType& parameters,
  ParameterType& momenta) const
{
  for(unsigned int i = 0; i < path_length; i++)
    LeapfrogStep(parameters, momenta);
}

template<class ParameterType>
void SimpleHamiltonian<ParameterType>::LeapfrogStep(
  ParameterType& parameters,
  ParameterType& momenta) const
{
  momenta -=
    this->model.EnergyPartials(parameters) * step_length / 2.;

  parameters +=
    momenta * step_length / getParameterMasses();

  momenta -=
    this->model.EnergyPartials(parameters) * step_length / 2.;
}

template<class ParameterType>
ParameterType SimpleHamiltonian<ParameterType>::RandomMomentum(const ParameterType& parameter, double temperature) const
{
  static thread_local std::random_device device;
  static thread_local std::mt19937_64 generator(device());
  static thread_local std::normal_distribution<double> normal(0., 1.);

  const double sqrt_temp = sqrt(temperature);

  ParameterType out;
  for (int i = 0; i < out.dimension; i++)
    out.parameters.at(i) = normal(generator) * 
      sqrt(getParameterMasses().parameters.at(i)) * sqrt_temp;

  return out;
}

#endif
