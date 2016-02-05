//// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
//// Released under the terms of the GNU General Public License version 3 or later.

//#include <gtest/gtest.h>

//#include "generic_step.hh"
//#include "mixins/eps.hh"
//#include "mixins/iterativeRefinements.hh"
//#include "mixins/verbosity.hh"

//namespace
//{
//  template <class Real>
//  struct MixinStep :
//      Dune::Mixin::IterativeRefinements ,
//      Dune::Mixin::Verbosity ,
//      Dune::Mixin::Eps<Real>
//  {};
//}

//TEST(GenericStep,Empty)
//{
//  auto step = Dune::GenericStep< MixinStep<double> >();
//  ASSERT_DOUBLE_EQ( step.getFinalIterate() , 0 );
//}

////TEST(GenericStep,Parameters)
////{
////  auto step = Dune::GenericStep<MixinStep>();
////  step.setIterativeRefinements(1);
////  ASSERT_EQ( step.getIterativeRefinements() , 1 );
////  step.setVerbosityLevel(2);
////  ASSERT_EQ( step.getVerbosityLevel() , 2 );
////  step.setEps( 1e-3 );
////  ASSERT_DOUBLE_EQ( step.eps() , 1e-3 );
////}
