#ifndef LENNARD_JONES_MODEL_H
#define LENNARD_JONES_MODEL_H

#include "parameter_set.inl"
#include "model.inl"
#include "point.inl"

template<int particles>
class LennardJonesParticles : public ParameterSet<particles * 2>
{
public:
  const Point getNthParticle(int n) const { return Point(&parameters.at(2 * n)); }
  void setNthParticle(int n, const Point& point)
  {
    parameters.at(2 * n) = point.at(0);
    parameters.at(2 * n + 1) = point.at(1);
  }

  static const int particle_count = particles;
  const int getParticleCount() const { return particles; }
};

template<int particles>
class LennardJonesModel : public Model<LennardJonesParticles<particles> >
{
public:
  LennardJonesModel(double particle_count, double x_size = 10, double y_size = 10);
private:
  double CalculateEnergy(const parameter_type& parameters) const = 0;
  parameter_type CalculateEnergyPartials(const parameter_type& parameters) const = 0;

  double DistanceEnergy(const Point& a, const Point& b) const;
  LennardJonesParticles<2> DistanceEnergyPartials(const Point& a, const Point& b) const;

  Point Displacement(const Point& a, const Point& b) const;

  double Distance(const Point& a, const Point& b) const;
  LennardJonesParticles<2> DistancePartials(const Point& a, const Point& b);
};

template <int particles>
LennardJonesModel<particles>::LennardJonesModel(
    double particle_count, double x_size, double y_size)
{

}

template <int particles>
double LennardJonesModel<particles>::CalculateEnergy(
    const LennardJonesModel<particles>::parameter_type& parameters) const
{
  double energy = 0.;
  for (int i = 0; i < parameter_type::particle_count; i++)
  {

  }
}

template <int particles>
LennardJonesModel<particles>::parameter_type 
LennardJonesModel<particles>::CalculateEnergyPartials(
    const LennardJonesModel<particles>::parameter_type& parameters) const
{

}

template <int particles>
double LennardJonesModel<particles>::DistanceEnergy(
    const Point& a, const Point& b)
{

}

template <int particles>
LennardJonesParticles<2> LennardJonesModel<particles>::DistanceEnergyPartials(
    const Point& a, const Point& b) const
{

}

template <int particles>
Point LennardJonesModel<particles>::Displacement(
    const Point& a, const Point& b)
{
  const Point half_bounds = 0.5 * bounds;
  Point difference = ( ( a - b ) % bounds + bounds ) % bounds;
  Point difference = ( ( difference + half_bounds ) % bounds ) - half_bounds;

  return difference;
}

template <int particles>
double LennardJonesModel<particles>::Distance(
    const Point& a, const Point& b)
{
  const Point displacement = Displacement(a, b);
  return displacement.Magnitude();
}

template <int particles>
LennardJonesParticles<2> LennardJonesModel<particles>::DistancePartials(
    const Point& a, const Point& b)
{

}

#endif
