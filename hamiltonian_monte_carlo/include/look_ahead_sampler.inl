#ifndef LOOK_AHEAD_SAMPLER_INL
#define LOOK_AHEAD_SAMPLER_INL

#include "hamiltonian_monte_carlo.inl"
#include <algorithm>

template<class ParameterType>
class LookAheadSampler<ParameterType> : public HamiltonianMonteCarlo<ParameterType>
{
private:
  int maximum_proposals = 50;

  class LookAheadProbability
  {
  private:
    int maximum_proposals;

    std::vector<std::vector<double> > dp_array;
    std::vector<double> point_probabilities;

    void setTransitionProb(int start, int end, double probability);

  public:
    LookAheadProbability(int maximum_proposals, double initial_prob)
      : maximum_proposals(maximum_proposals) 
    {
      dp_array.reserve(maximum_proposals + 1);

      point_probabilities.reserve(maximum_proposals + 1);
      point_probabilities.push_back(initial_prob);
    }

    double getTransitionCount()
    {
      return point_probabilities - 1;
    }

    double getTransitionProb(int start, int end);

    double getNextProbability(double point_prob)
    {
      point_probabilities.push_back(point_prob);
      return getTransitionProb(0, getTransitionCount());
    }
 
    double getNthTransitionProbability(int transition)
    {
      return getTransitionProb(0, transition + 1);
    }
  };

public:
  LookAheadSampler(const Model& model, const Hamiltonian& hamiltonian)
    : HamiltonianMonteCarlo(model, hamiltonian) {}

  virtual void SimulateStep(ParameterType& parameter) override;

  void setMaximumProposals(int new_maximum) { maximum_proposals = new_maximum; }
  int getMaximumProposals() { return maximum_proposals; }
};

template<class ParameterType>
void LookAheadSampler<ParameterType>::LookAheadProbability::setTransitionProb(int start, int end, double probability)
{
  if (dp_array.size() < start)
  {
    const std::vector<double> row;
    row.reserve(end);
    for (;dp_array.size() < start;)
      dp_array.push_back(row)
  }

  if (dp_array.at(start).size() < end)
  {
    for (;dp_array.at(start).size() < end)
      dp_array.at(start).push_back(-1.);
  }

  dp_array.at(start).at(end) = probability;
}

template<class ParameterType>
double LookAheadSampler<ParameterType>::LookAheadProbability::getTransitionProb(int start, int end)
{
  if (!(dp_array.size() > start &&
      dp_array.at(start).size() > end &&
      dp_array.at(start).at(end) >= 0.))
  {
    int direction = (start < end) - (end < start);

    double forwards_remaining_prob = 1.;
    for (int i = start + direction; i != end; i += direction)
    {
      forwards_remaining_prob -= getTransitionProb(start, i);
    }

    double backwards_remaining_prob = 1.;
    for (int i = end - direction; i != start; i -= direction)
    {
      backwards_remaining_prob -= getTransitionProb(end, i);
    }
    backwards_remaining_prob *=
      point_probabilities.at(end) / point_probabilities.at(start);

    const double calculated_prob = std::max(0, 
      std::min(backwards_remaining_prob, forwards_remaining_prob));

    setTransitionProb(int start, int end, calculated_prob);
  }
  return dp_array.at(start).at(end);
}

template<class ParameterType>
void LookAheadSampler<ParameterType>::SimulateStep(ParameterType& parameter)
{
  //Choose the cutoff probability
  double cutoff = MarkovChainMonteCarlo::getRandomUniform();
  double cumulative_probability = 0.;

  //Create dynamic programming object
  LookAheadProbability prob_calc(
    maximum_proposals,
    exp(model.Energy(parameter)));

  ParameterType old_parameter = parameter;

  for (int i = 0; i < maximum_proposals; i++)
  {
    ParameterType momentum = hamiltonian.RandomMomentum();
    hamiltonian.GenerateStep(parameter, momentum);
  
    cumulative_probability += 
      getNextProbability(exp(model.Energy(parameter)));

    if (cumulative_probability > cutoff)
      return;
  }

  parameter = old_parameter;
  return;
}

#endif
