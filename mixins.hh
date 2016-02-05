#ifndef DUNE_MIXINS_HH
#define DUNE_MIXINS_HH

#include "mixins/absoluteAccuracy.hh"
#include "mixins/eps.hh"
#include "mixins/iterativeRefinements.hh"
#include "mixins/maxSteps.hh"
#include "mixins/minimalAccuracy.hh"
#include "mixins/relativeAccuracy.hh"
#include "mixins/verbosity.hh"

#define DUNE_ISTL_MIXINS(Real) Mixin::AbsoluteAccuracy<Real>, Mixin::MinimalAccuracy<Real>, Mixin::RelativeAccuracy<Real>, Mixin::Verbosity, \
  Mixin::Eps<Real>, Mixin::IterativeRefinements, Mixin::MaxSteps

#endif // DUNE_MIXINS_HH

