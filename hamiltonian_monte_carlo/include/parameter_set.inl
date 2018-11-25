#ifndef PARAMETER_SET_INL
#define PARAMETER_SET_INL

#include <array>

template<int N>
class ParameterSet<N>
{
public:
  static const int dimension = N;
  std::array<double, N> parameters;

  explicit ParameterSet<N>(double initial_values = 0.)
  {
    for (int i = 0; i < N; i++)
      parameters.at(i) = initial_values;
  }
  explicit ParameterSet<N>(std::array<double, N> initial_values)
  {
    parameters = initial_values;
  }

  ParameterSet<N> operator+(const ParameterSet<N> other) const
  {
    ParameterSet<N> out = *this;
    out += other;
    return out;
  }

  ParameterSet<N> operator*(const ParameterSet<N> other) const
  {
    ParameterSet<N> out = *this;
    out *= other;
    return out;
  }

  ParameterSet<N> operator-(const ParameterSet<N> other)
  {
    ParameterSet<N> out = *this;
    out -= other;
    return out;
  }

  ParameterSet<N> operator/(const ParameterSet<N> other)
  {
    ParameterSet<N> out = *this;
    out /= other;
    return out;
  }

  ParameterSet<N>& operator+=(const ParameterSet<N> other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) += other.at(i);

    return *this;
  }

  ParameterSet<N>& operator*=(const ParameterSet<N> other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) *= other.at(i);

    return *this;
  }

  ParameterSet<N>& operator-=(const ParameterSet<N> other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) -= other.at(i);

    return *this;
  }

  ParameterSet<N>& operator/=(const ParameterSet<N> other)
  {
    for(int i = 0; i < N; i++)
      parameters.at(i) /= other.at(i);

    return *this;
  }

};

#endif
