#ifndef TRAJECTORY_TEMPERED_HAMILTONIAN_INL
#define TRAJECTORY_TEMPERED_HAMILTONIAN_INL

#include "simple_hamiltonian.inl"

template<class ParameterType>
class TrajectoryTemperedHamiltonian<ParameterType> : public SimpleHamiltonian<ParameterType>
{
private:
  double alpha_ratio;

public:
  explicit TrajectoryTemperedHamiltonian(
    const Model& model,
    double step_length = 0.01,
    int path_length = 100,
    double alpha_ratio = 1.05)
    : SimpleHamiltonian(model, step_length, path_length)
      alpha_ratio(alpha_ratio) {}

  double getAlphaRatio() const { return alpha_ratio; }
  void setAlphaRatio(double new_ratio) { alpha_ratio = new_ratio; }

  virtual void IntegratePath(
    ParameterType parameters&,
    ParameterType momenta&) const override;
};

template<class ParameterType>
void TrajectoryTemperedHamiltonian<ParameterType>::IntegratePath(
  ParameterType parameters&,
  ParameterType momenta&) const
{
  const double halfway = path_length / 2.;
  for(unsigned int i = 0; i < path_length; i++)
  {
    if (i < halfway)
      momenta *= alpha_ratio;
    else if (i > halfway)
      momenta /= alpha_ratio;
    LeapfrogStep(parameters, momenta);
  }
}

#endif
