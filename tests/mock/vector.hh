#ifndef DUNE_ISTL_TESTS_MOCK_VECTOR_HH
#define DUNE_ISTL_TESTS_MOCK_VECTOR_HH

#include <dune/common/typetraits.hh>

namespace Dune
{
  namespace Mock
  {
    struct Vector
    {
      using field_type = double;

      Vector();

      Vector& operator+=(const Vector&);
    };
  }

  template <>
  struct FieldTraits<Mock::Vector>
  {
    using field_type = double;
    using real_type = double;
  };
}

#endif // DUNE_ISTL_TESTS_MOCK_VECTOR_HH
