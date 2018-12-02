#ifndef LENNARD_JONES_MODEL_INL
#define LENNARD_JONES_MODEL_INL

#include "parameter_set.inl"
#include "model.inl"
#include "point.inl"
#include "poisson_generator.inl"

template<int particles>
class LennardJonesParticles : public ParameterSet<particles * 2>
{
public:
  LennardJonesParticles() {}
  LennardJonesParticles(double in) : ParameterSet<particles * 2>(in) {}
  LennardJonesParticles(const ParameterSet<particles * 2>& old) 
    : ParameterSet<particles * 2>(old) {}

  Point getNthParticle(int n) const { return Point(&this->parameters.at(2 * n)); }
  void setNthParticle(int n, const Point& point)
  {
    this->parameters.at(2 * n) = point.parameters.at(0);
    this->parameters.at(2 * n + 1) = point.parameters.at(1);
  }

  static const int particle_count = particles;
  const int getParticleCount() const { return particles; }
};

template<int particles>
class LennardJonesModel : public Model<LennardJonesParticles<particles> >
{
private:
  Point bounds;

  double potential_depth = 1.;
  double zero_potential_distance = 1.;
  double truncation_distance = 3.;

  double continuity_adjustment_value;

public:
  LennardJonesModel(double x_size = 10, double y_size = 10);

  using parameter_type = typename Model<LennardJonesParticles<particles> >::parameter_type;

  virtual parameter_type ParameterMapReals(const parameter_type& parameter) const override;

  parameter_type getRandomInitialState() const override;

  static std::vector<double> CalculateMeanSquaredDisplacement(
      const std::vector<parameter_type>& in);
private:
  double CalculateEnergy(const parameter_type& parameters) const;
  parameter_type CalculateEnergyPartials(const parameter_type& parameters) const;

  double DistanceEnergy(const Point& a, const Point& b) const;
  LennardJonesParticles<2> DistanceEnergyPartials(const Point& a, const Point& b) const;

  Point Displacement(const Point& a, const Point& b) const;

  double Distance(const Point& a, const Point& b) const;
  LennardJonesParticles<2> DistancePartials(const Point& a, const Point& b) const;

  double Potential(double distance) const;
  double PotentialPartials(double distance) const;

  virtual parameter_type RealMapPartials(
      const parameter_type& mapped_parameters,
      const parameter_type& partials) const override
  { return partials; }
};

template <int particles>
LennardJonesModel<particles>::LennardJonesModel(
    double x_size, double y_size)
{
  bounds.setX(x_size);
  bounds.setY(y_size);

  continuity_adjustment_value = 0;
  continuity_adjustment_value =
    Potential(truncation_distance);
}

template <int particles>
typename LennardJonesModel<particles>::parameter_type 
LennardJonesModel<particles>::ParameterMapReals(
    const LennardJonesModel<particles>::parameter_type& parameter) const
{
  parameter_type out;
  for (int i = 0; i < particles; i++)
    out.setNthParticle(
        i,
        ( ( out.getNthParticle(i) % bounds ) + bounds ) % bounds);
  return out;
}

template <int particles>
double LennardJonesModel<particles>::CalculateEnergy(
    const LennardJonesModel<particles>::parameter_type& parameters) const
{
  double energy = 0.;
  for (int i = 0; i < parameter_type::particle_count; i++)
  {
    for (int j = i + 1; j < parameter_type::particle_count; j++)
    {
      energy +=
        DistanceEnergy(
            parameters.getNthParticle(i),
            parameters.getNthParticle(j));
    }
  }
  return energy;
}

template <int particles>
typename LennardJonesModel<particles>::parameter_type 
LennardJonesModel<particles>::CalculateEnergyPartials(
    const LennardJonesModel<particles>::parameter_type& parameters) const
{
  parameter_type partials;
  for (int i = 0; i < parameter_type::particle_count; i++)
  {
    for (int j = i + 1; j < parameter_type::particle_count; j++)
    {
      LennardJonesParticles<2> tmp_partials = 
        DistanceEnergyPartials(
            parameters.getNthParticle(i),
            parameters.getNthParticle(j));
      partials.setNthParticle(
          i,
          partials.getNthParticle(i) + tmp_partials.getNthParticle(0) );
      partials.setNthParticle(
          j,
          partials.getNthParticle(j) + tmp_partials.getNthParticle(1) );
    }
  }
  return partials;
}

template <int particles>
double LennardJonesModel<particles>::DistanceEnergy(
    const Point& a, const Point& b) const
{
  const double distance = Distance(a, b);
  return Potential(distance);
}

template <int particles>
LennardJonesParticles<2> LennardJonesModel<particles>::DistanceEnergyPartials(
    const Point& a, const Point& b) const
{
  const double distance = Distance(a, b);
  const LennardJonesParticles<2> dist_partials =
    DistancePartials(a, b);
  const double potential_partial =
    PotentialPartials(distance);

  const LennardJonesParticles<2> partials =
    dist_partials * potential_partial;
  return partials;
}

template <int particles>
Point LennardJonesModel<particles>::Displacement(
    const Point& a, const Point& b) const
{
  const Point half_bounds = bounds * 0.5;
  Point difference = ( ( a - b ) + bounds ) % bounds;
  difference = ( ( difference + half_bounds ) % bounds ) - half_bounds;

  return difference;
}

template <int particles>
double LennardJonesModel<particles>::Distance(
    const Point& a, const Point& b) const
{
  const Point displacement = Displacement(a, b);
  return displacement.Magnitude();
}

template <int particles>
LennardJonesParticles<2> LennardJonesModel<particles>::DistancePartials(
    const Point& a, const Point& b) const
{
  LennardJonesParticles<2> partials;
  const Point difference = Displacement(a, b);
  const double distance = difference.Magnitude();
  const Point ratio = difference / distance;
  partials.setNthParticle(0, ratio);
  partials.setNthParticle(1, ratio * -1);
  return partials;
}

template <int particles>
double LennardJonesModel<particles>::Potential(double distance) const
{
  if (distance > truncation_distance)
    return 0.;

  const double recip_distance = zero_potential_distance / distance;
  const double recip_squared = recip_distance * recip_distance;
  const double recip_six = 
    recip_squared * recip_squared * recip_squared;
  const double recip_twelve = recip_six * recip_six;

  const double potential =
    4 * potential_depth * (recip_twelve - recip_six);

  const double final_potential = 
    potential - continuity_adjustment_value;

  return potential;
}

template <int particles>
double LennardJonesModel<particles>::PotentialPartials(double distance) const
{
  if (distance > truncation_distance)
    return 0.;

  const double recip_distance = zero_potential_distance / distance;
  const double recip_squared = recip_distance * recip_distance;
  const double recip_six = 
    recip_squared * recip_squared * recip_squared;
  const double recip_twelve = recip_six * recip_six;

  const double partial =
    24. * potential_depth * recip_distance *
    (recip_six - 2. * recip_twelve);
  return partial;
}

template<int particles>
typename LennardJonesModel<particles>::parameter_type
LennardJonesModel<particles>::getRandomInitialState() const
{
  PoissonGenerator::DefaultPRNG generator;

  parameter_type out;

  //number of points, prng, test point count, fill circle/rectangle, minimum distance
  std::vector<PoissonGenerator::sPoint> points = 
        PoissonGenerator::GeneratePoissonPoints(
      particles,
      generator,
      30,
      false,
      0.5);

  int index = 0;
  for (PoissonGenerator::sPoint point : points)
  {
    Point tmp_point = Point(point.x, point.y);
    out.setNthParticle(index, tmp_point);
    index++;
  }

  return out;
}

template<int particles>
std::vector<double> LennardJonesModel<particles>::CalculateMeanSquaredDisplacement(
    const std::vector<typename LennardJonesModel<particles>::parameter_type>& in)
{
  std::vector<double> out;

  for (int i = 0; i < in.size(); i++)
  {
    double total_squared_displacement;
    for (int j = 0; j < particles; j++)
    {
      total_squared_displacement +=
        (in.at(0).getNthParticle(j) -
         in.at(i).getNthParticle(j)).SquaredMagnitude();
    }
    out.push_back(total_squared_displacement / particles);
  }
  return out;
}
#endif
