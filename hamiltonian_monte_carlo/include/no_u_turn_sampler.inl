#ifndef NO_U_TURN_SAMPLER_INL
#define NO_U_TURN_SAMPLER_INL

#include "hamiltonian_monte_carlo.inl"

template<class ParameterType>
class NUTS<ParameterType> : public HamiltonianMonteCarlo
{

};

#endif
