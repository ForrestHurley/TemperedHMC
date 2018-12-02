#ifndef NO_U_TURN_SAMPLER_INL
#define NO_U_TURN_SAMPLER_INL

#include "hamiltonian_monte_carlo.inl"
#include "simple_hamiltonian.inl"
#include <cassert>

template<class ParameterType>
class NUTS : public HamiltonianMonteCarlo<ParameterType>
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
  NUTS(const Model<ParameterType>& model, const SimpleHamiltonian<ParameterType>& hamiltonian)
    : HamiltonianMonteCarlo<ParameterType>(model, hamiltonian)
  {
    assert(hamiltonian.getPathLength() == 1);
  }

  virtual void SimulateStep(ParameterType& parameter) override;

  void setMaximumDeltaEnergy(double new_max) { maximum_delta_energy = new_max; }
  double getMaximumDeltaEnergy() const { return maximum_delta_energy; }

};

template<class ParameterType>
typename NUTS<ParameterType>::TreeReturn NUTS<ParameterType>::BuildTree(
    const ParameterType& parameter,
    const ParameterType& momentum,
    double cutoff,
    int direction,
    int height)
{
  if (height == 0)
  {
    ParameterType tmp_parameter = parameter;
    ParameterType tmp_momentum = momentum;
    //Base case
    tmp_momentum *= direction;
    this->hamiltonian.GenerateStep(tmp_parameter, tmp_momentum);
    tmp_momentum *= direction;

    TreeReturn out;
    out.left_parameter = tmp_parameter;
    out.left_momentum = tmp_momentum;
    out.right_parameter = tmp_parameter;
    out.right_momentum = tmp_momentum;
    out.new_parameter = tmp_parameter;

    const double energy = this->hamiltonian.Energy(tmp_parameter, tmp_momentum);

    out.count = (cutoff <= exp(-energy));
    out.safe = (cutoff < exp(-energy + maximum_delta_energy));

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
    }
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
    if (this->getRandomUniform() < second_param_prob)
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
  ParameterType momentum = this->hamiltonian.RandomMomentum(parameter);

  const double initial_energy = this->hamiltonian.Energy(parameter, momentum);
  double initial_probability = exp(-initial_energy);
  double cutoff = this->getRandomUniform() * initial_probability;  

  TreeReturn working_tree;
  working_tree.left_parameter = parameter;
  working_tree.left_momentum = momentum;
  working_tree.right_parameter = parameter;
  working_tree.right_momentum = momentum;
  working_tree.new_parameter = parameter;
  working_tree.count = 1;
  working_tree.safe = true;

  bool made_step = false;

  unsigned int max_iters = 50000;
  unsigned int height = 0;
  while (working_tree.safe == true & height < max_iters)
  {

    int direction = 
      (this->getRandomUniform() > 0.5) * 2 - 1;

    TreeReturn second_tree;
    if (direction == -1)
    {
      second_tree = BuildTree(
        working_tree.left_parameter, 
        working_tree.left_momentum, 
        cutoff, direction, height);
  
      working_tree.left_parameter = second_tree.left_parameter;
      working_tree.left_momentum = second_tree.left_momentum;      
    }
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
        (this->getRandomUniform() 
         < (double)second_tree.count / working_tree.count))
    {
      parameter = second_tree.new_parameter;
      made_step = true;
    }

    const ParameterType difference = working_tree.right_parameter - working_tree.left_parameter;
    bool expanding_right =  
      (difference * working_tree.right_momentum).Sum()
      >= 0;
    bool expanding_left =  
      (difference * working_tree.left_momentum).Sum()
      >= 0;

    working_tree.safe =
      second_tree.safe & expanding_right & expanding_left;
    working_tree.count += second_tree.count;

    height++;
  }

  if (made_step)
    this->accept_count++;
  else
    this->reject_count++;
}

#endif
