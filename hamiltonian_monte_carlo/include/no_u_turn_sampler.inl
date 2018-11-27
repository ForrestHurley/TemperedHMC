#ifndef NO_U_TURN_SAMPLER_INL
#define NO_U_TURN_SAMPLER_INL

#include "hamiltonian_monte_carlo.inl"
#include "simple_hamiltonian.inl"
#include <cassert>

template<class ParameterType>
class NUTS<ParameterType> : public HamiltonianMonteCarlo<ParameterType>
{
private:
  struct TreeReturn {
    ParameterType left_parameter;
    ParameterType left_momentum;
    ParameterType right_parameter;
    ParameterType right_momentum;
    ParameterType new_parameter;
    int count;
    bool safe;
  };

  TreeReturn BuildTree(
    const ParameterType& parameter,
    const ParameterType& momentum,
    double cutoff,
    int direction,
    int height);

  double maximum_delta_energy = 1000;

public:
  NUTS(const Model& model, const SimpleHamiltonian& hamiltonian)
    : HamiltonianMonteCarlo(model, hamiltonian)
  {
    assert(hamiltonian.getPathLength() == 1);
  }

  virtual void SimulateStep(ParameterType& parameter) override;

  void setMaximumDeltaEnergy(double new_max) { maximum_delta_energy = new_max; }
  double getMaximumDeltaEnergy() const { return maximum_delta_energy; }

};

template<class ParameterType>
NUTS<ParameterType>::TreeReturn NUTS<ParameterType>::BuildTree(
    const ParameterType& parameter,
    const ParameterType& momentum,
    double cutoff,
    int direction,
    int height)
{
  if (height == 0)
  {
    //Base case
    momentum *= direction;
    hamiltonian.GenerateStep(parameter, momentum);
    momentum *= direction;

    TreeReturn out;
    out.left_parameter = parameter;
    out.left_momentum = momentum;
    out.right_parameter = parameter;
    out.right_momentum = momentum;
    out.new_parameter = parameter;

    double energy = model.Energy(parameter) + 
        static_cast<const SimpleHamiltonian&>(hamiltonian).KineticEnergy(momentum);

    out.count = (u <= exp(-energy));
    out.safe = (u < exp(-energy + maximum_delta_energy));

    return out;
  }

  TreeReturn first_tree =
    BuildTree(parameter, momentum, cutoff, direction, height - 1);

  TreeReturn second_tree;
  if (first_tree.safe)
  {
    if (direction == -1)
    {
      second_tree = BuildTree(
        first_tree.left_parameter, 
        first_tree.left_momentum, 
        cutoff, direction, height - 1);
  
      first_tree.left_parameter = second_tree.left_parameter;
      first_tree.left_momentum = second_tree.left_momentum;      
    {
    if (direction == 1)
    {
      second_tree = BuildTree(
        first_tree.right_parameter,
        first_tree.right_momentum,
        cutoff, direction, height - 1);

      first_tree.right_parameter = second_tree.right_parameter;
      first_tree.right_momentum = second_tree.right_momentum;
    }

    double second_param_prob = 
      (double) second_tree.count / ( first_tree.count + second_tree.count);
    if (getRandomUniform() < second_param_prob)
      first_tree.new_parameter = second_tree.new_parameter;

    const ParameterType difference = first_tree.right_parameter - first_tree.left_parameter;
    bool expanding_right =  
      (difference * first_tree.right_momentum).Sum()
      >= 0;
    bool expanding_left =  
      (difference * first_tree.left_momentum).Sum()
      >= 0;

    first_tree.safe =
      second_tree.safe & expanding_right & expanding_left;
    first_tree.count += second_tree.count;
  }

  return first_tree;
}

template<class ParameterType>
void NUTS<ParameterType>::SimulateStep(ParameterType& parameter)
{
  ParameterType momentum = hamiltonian.RandomMomentum();

  double initial_energy = model.Energy(parameter) +
    static_cast<const SimpleHamiltonian&>(hamiltonian).KineticEnergy(momentum);
  double initial_probability = exp(-initial_energy);
  double u = MarkovChainMonteCarlo::getRandomUniform() * intial_probability;  

  TreeReturn working_tree;
  working_tree.left_parameter = parameter;
  working_tree.left_momentum = momentum;
  working_tree.right_parameter = parameter;
  working_tree.right_momentum = momentum;
  working_tree.new_parameter = parameter;
  working_tree.count = 1;
  working_tree.safe = true;

  unsigned int max_iters = 50000;
  unsigned int j = 0;
  while (working_tree.safe == true & j < max_iters)
  {

    int direction = 
      (MarkovChainMonteCarlo::getRandomUniform() > 0.5) * 2 - 1;

    TreeReturn second_tree;
    if (direction == -1)
    {
      second_tree = BuildTree(
        working_tree.left_parameter, 
        working_tree.left_momentum, 
        cutoff, direction, height);
  
      working_tree.left_parameter = second_tree.left_parameter;
      working_tree.left_momentum = second_tree.left_momentum;      
    {
    if (direction == 1)
    {
      second_tree = BuildTree(
        working_tree.right_parameter,
        working_tree.right_momentum,
        cutoff, direction, height);

      working_tree.right_parameter = second_tree.right_parameter;
      working_tree.right_momentum = second_tree.right_momentum;
    }

    if (second_tree.safe && 
        (MarkovChainMonteCarlo::getRandomUniform() 
         < (double)second_tree.count / working_tree.count))
    {
      parameter = second_tree.new_parameter;
    }

    const ParameterType difference = first_tree.right_parameter - first_tree.left_parameter;
    bool expanding_right =  
      (difference * first_tree.right_momentum).Sum()
      >= 0;
    bool expanding_left =  
      (difference * first_tree.left_momentum).Sum()
      >= 0;

    working_tree.safe =
      second_tree.safe & expanding_right & expanding_left;
    working_tree.count += second_tree.count;

    j++;
  }
}

#endif
