#ifndef PARAMETER_SET_INL
#define PARAMETER_SET_INL

#include <array>

template<int N>
class ParameterSet<N>
{
public:
  static const int dimension = N;
  std::array<double, N> parameters;

  ParameterSet<N>()
  {
    for (int i = 0; i < N; i++)
      parameters.at(i) = 0.;
  }
  explicit ParameterSet<N>(double initial_values = 0.)
  {
    for (int i = 0; i < N; i++)
      parameters.at(i) = initial_values;
  }
  explicit ParameterSet<N>(std::array<double, N> initial_values)
  {
    parameters = initial_values;
  }

  ParameterSet<N>& operator%=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) %= other.at(i);

    return *this;
  }

  ParameterSet<N>& operator+=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) += other.at(i);

    return *this;
  }

  ParameterSet<N>& operator*=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) *= other.at(i);

    return *this;
  }

  ParameterSet<N>& operator-=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) -= other.at(i);

    return *this;
  }

  ParameterSet<N>& operator/=(const ParameterSet<N>& other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) /= other.at(i);

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

  ParameterSet<N> operator-(const ParameterSet<N>& other)
  {
    ParameterSet<N> out = *this;
    out -= other;
    return out;
  }

  ParameterSet<N> operator/(const ParameterSet<N>& other)
  {
    ParameterSet<N> out = *this;
    out /= other;
    return out;
  }

  PrameterSet<N> operator%(const ParameterSet<N>& other)
  {
    ParameterSet<N> out = *this;
    out %= other;
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
    double sum = 0.
    for (double value : parameters)
      sum += value * value;
    return sum;
  }
  
  double Magnitude() const
  {
    return sqrt(MagnitudeSquared());
  }

};

#endif
