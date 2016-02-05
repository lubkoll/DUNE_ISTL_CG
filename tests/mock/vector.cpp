#include "vector.hh"

namespace Dune
{
  namespace Mock
  {
    Vector::Vector()
    {}

    Vector& Vector::operator+=(const Vector&)
    {
      return *this;
    }
  }
}
