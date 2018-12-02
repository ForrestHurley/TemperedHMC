#ifndef PARAMETER_SET_INL
#define PARAMETER_SET_INL

#include <array>
#include "math.h"

//Do not inherit from this
//Inherit from ParameterSet
class ParameterBase {};

template<int N>
class ParameterSet : ParameterBase
{
public:
  static const int dimension = N;
  std::array<double, N> parameters;

  ParameterSet<N>()
  {
    for (int i = 0; i < N; i++)
      parameters.at(i) = 0.;
  }
  explicit ParameterSet<N>(double initial_values)
  {
    for (int i = 0; i < N; i++)
      parameters.at(i) = initial_values;
  }
  explicit ParameterSet<N>(std::array<double, N> initial_values)
  {
    parameters = initial_values;
  }
  ParameterSet<N>(const ParameterSet<N>& old)
  {
    parameters = old.parameters;
  }

  ParameterSet<N>& operator%=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) = std::fmod(
          parameters.at(i),
          other.parameters.at(i));

    return *this;
  }

  ParameterSet<N>& operator+=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) += other.parameters.at(i);

    return *this;
  }

  ParameterSet<N>& operator*=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) *= other.parameters.at(i);

    return *this;
  }

  ParameterSet<N>& operator-=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) -= other.parameters.at(i);

    return *this;
  }

  ParameterSet<N>& operator/=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) /= other.parameters.at(i);

    return *this;
  }

  ParameterSet<N>& operator+=(double other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) += other;

    return *this;
  }

  ParameterSet<N>& operator*=(double other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) *= other;

    return *this;
  }

  ParameterSet<N>& operator-=(double other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) -= other;

    return *this;
  }


  ParameterSet<N>& operator/=(double other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) /= other;

    return *this;
  }

  ParameterSet<N> operator+(const ParameterSet<N>& other) const
  {
    ParameterSet<N> out = *this;
    out += other;
    return out;
  }

  ParameterSet<N> operator*(const ParameterSet<N>& other) const
  {
    ParameterSet<N> out = *this;
    out *= other;
    return out;
  }

  ParameterSet<N> operator-(const ParameterSet<N>& other) const
  {
    ParameterSet<N> out = *this;
    out -= other;
    return out;
  }

  ParameterSet<N> operator/(const ParameterSet<N>& other) const
  {
    ParameterSet<N> out = *this;
    out /= other;
    return out;
  }

  ParameterSet<N> operator%(const ParameterSet<N>& other) const
  {
    ParameterSet<N> out = *this;
    out %= other;
    return out;
  }

  ParameterSet<N> operator+(double other) const
  {
    ParameterSet<N> out = *this;
    out += other;
    return out;
  }

  ParameterSet<N> operator*(double other) const
  {
    ParameterSet<N> out = *this;
    out *= other;
    return out;
  }

  ParameterSet<N> operator-(double other) const
  {
    ParameterSet<N> out = *this;
    out -= other;
    return out;
  }

  ParameterSet<N> operator/(double other) const
  {
    ParameterSet<N> out = *this;
    out /= other;
    return out;
  }

  double Sum() const
  {
    double sum = 0.;
    for (double value : parameters)
      sum += value;
    return sum;
  }

  double MagnitudeSquared() const
  {
    double sum = 0.;
    for (double value : parameters)
      sum += value * value;
    return sum;
  }
  
  double Magnitude() const
  {
    return sqrt(MagnitudeSquared());
  }

  ParameterSet<N> Sign() const
  {
    ParameterSet<N> out;
    for (int i = 0; i < N; i++)
      out.at(i) = 
        (parameters.at(i) > 0) -
        (parameters.at(i) < 0);
    return out;
  }
};

#endif
