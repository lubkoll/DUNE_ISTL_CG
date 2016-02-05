// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#include "cg_solver.hh"

struct Vector
{
  double value = 0;
};

struct Operator : public Dune::LinearOperator<Vector,Vector>
{
  void apply(const Vector& x, Vector& y) const final override
  {
    y = x;
  }
};


