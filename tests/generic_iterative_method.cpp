// Copyright (C) 2015 by Lars Lubkoll. All rights reserved.
// Released under the terms of the GNU General Public License version 3 or later.

#include <string>

#include <gtest/gtest.h>

#include <dune/istl/solvers.hh>

#include "../generic_iterative_method.hh"

#include "mock/step.hh"
#include "mock/termination_criteria.hh"
#include "mock/vector.hh"

namespace Mock = Dune::Mock;

using Mock::Step;
using Mock::RestartingStep;
using Mock::TerminatingStep;
using Mock::TerminationCriterion;
using Mock::MixinTerminationCriterion;

namespace
{
  inline double testAccuracy()
  {
    return 1e-6;
  }
}


TEST(GenericIterativeMethod,NotConvergedInZeroSteps)
{
  auto dummySolver = Dune::makeGenericIterativeMethod(Step(),TerminationCriterion<Step>());
  dummySolver.setMaxSteps(0);
  Mock::Vector x, b;

  EXPECT_FALSE( dummySolver.wasInitialized );
  Dune::InverseOperatorResult info;
  dummySolver.apply(x,b,info);
  EXPECT_FALSE( info.converged );
  EXPECT_TRUE( dummySolver.wasInitialized );
  EXPECT_TRUE( dummySolver.getTerminationCriterion().wasInitialized );
  EXPECT_FALSE( dummySolver.wasReset );
}

TEST(GenericIterativeMethod,NotConverged)
{
  auto dummySolver = Dune::makeGenericIterativeMethod(Step(),TerminationCriterion<Step>(false));
  dummySolver.setMaxSteps(10);
  Mock::Vector x, b;

  EXPECT_FALSE( dummySolver.wasInitialized );
  Dune::InverseOperatorResult info;
  dummySolver.apply(x,b,info);
  EXPECT_FALSE( info.converged );
  EXPECT_TRUE( dummySolver.wasInitialized );
  EXPECT_TRUE( dummySolver.getTerminationCriterion().wasInitialized );
  EXPECT_FALSE( dummySolver.wasReset );
}

TEST(GenericIterativeMethod,Converged)
{
  auto dummySolver = Dune::makeGenericIterativeMethod(Step(),TerminationCriterion<Step>());
  dummySolver.setMaxSteps(2);
  Mock::Vector x, b;

  EXPECT_FALSE( dummySolver.wasInitialized );
  Dune::InverseOperatorResult info;
  dummySolver.apply(x,b,info);
  EXPECT_TRUE( info.converged );
  EXPECT_TRUE( dummySolver.wasInitialized );
  EXPECT_TRUE( dummySolver.getTerminationCriterion().wasInitialized );
  EXPECT_FALSE( dummySolver.wasReset );
}

TEST(GenericIterativeMethod,Converged_TerminatingStep)
{
  TerminatingStep step;
  step.doTerminate = true;
  auto dummySolver = Dune::makeGenericIterativeMethod(step,TerminationCriterion<TerminatingStep>());
  dummySolver.setMaxSteps(1);
  Mock::Vector x, b;

  Dune::InverseOperatorResult info;
  dummySolver.apply(x,b,info);
  EXPECT_TRUE( info.converged );
  EXPECT_TRUE( dummySolver.wasInitialized );
  EXPECT_FALSE( dummySolver.wasReset );
}

TEST(GenericIterativeMethod,RestartAndTerminate)
{
  RestartingStep step(true);
  auto dummySolver = Dune::makeGenericIterativeMethod(step,TerminationCriterion<RestartingStep>(false));
  dummySolver.setMaxSteps(1);
  Mock::Vector x, b;

  Dune::InverseOperatorResult info;
  dummySolver.apply(x,b,info);
  EXPECT_TRUE( info.converged );
  EXPECT_TRUE( dummySolver.wasInitialized );
  EXPECT_TRUE( dummySolver.wasReset );
}

TEST(GenericIterativeMethod,MixinParameters)
{
  auto dummySolver = Dune::makeGenericIterativeMethod(Step(),MixinTerminationCriterion<Step>());
  dummySolver.setAbsoluteAccuracy(testAccuracy());
  EXPECT_DOUBLE_EQ( dummySolver.getTerminationCriterion().absoluteAccuracy() , testAccuracy() );
  dummySolver.setRelativeAccuracy(testAccuracy());
  EXPECT_DOUBLE_EQ( dummySolver.getTerminationCriterion().relativeAccuracy() , testAccuracy() );
  dummySolver.setMinimalAccuracy(testAccuracy());
  EXPECT_DOUBLE_EQ( dummySolver.getTerminationCriterion().minimalAccuracy() , testAccuracy() );
  dummySolver.setEps(testAccuracy());
  EXPECT_DOUBLE_EQ( dummySolver.getTerminationCriterion().eps() , testAccuracy() );
  dummySolver.setVerbosityLevel(2);
  EXPECT_EQ( dummySolver.getTerminationCriterion().verbosityLevel() , 2 );
}
