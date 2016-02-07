#ifndef DUNE_ISTL_TESTS_MOCK_LINEAR_OPERATOR_2D_H
#define DUNE_ISTL_TESTS_MOCK_LINEAR_OPERATOR_2D_H

#include <vector>
#include <dune/istl/operators.hh>

#include "vector.hh"

namespace Dune
{
  namespace Mock
  {
    class LinearOperator_2d : public LinearOperator<Vector,Vector>
    {
    public:
      LinearOperator_2d();

      void apply( const Vector& x, Vector& y ) const;

      void applyscaleadd( double a, const Vector& x, Vector& y ) const;

    private:
      std::vector< std::vector< double > > data_;
    };
  }
}

#endif // DUNE_ISTL_TESTS_MOCK_LINEAR_OPERATOR_2D_H
