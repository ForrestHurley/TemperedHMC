#ifndef POINT_INL
#define POINT_INL

#include "parameter_set.inl"

class Point : public ParameterSet<2>
{
public:
  Point() {}
  Point(double in)
  {
    parameters.at(0) = in;
    parameters.at(1) = in;
  }
  Point(double x, double y)
  {
    parameters.at(0) = x;
    parameters.at(1) = y;
  }
  Point(const double* array)
  {
    parameters.at(0) = *array;
    parameters.at(1) = *(array + 1);
  }
  Point(const ParameterSet<2>& old) : ParameterSet<2>(old) {}

  double& getX() { return parameters.at(0); }
  double getX() const { return parameters.at(0); }
  void setX(double new_x) { parameters.at(0) = new_x; }

  double& getY() { return parameters.at(1); }
  double getY() const { return parameters.at(1); }
  void setY(double new_y) { parameters.at(1) = new_y; }

};

#endif
