#include <gtest/gtest.h>

#include "dune/istl/solvers.hh"
#include "../relative_energy_termination_criterion.hh"

#include "mock/step.hh"

namespace
{
  struct TestRelativeEnergyErrorCriterion : ::testing::Test
  {
    TestRelativeEnergyErrorCriterion()
    {
      terminationCriterion.connect(step);
      terminationCriterion.init();
    }

    Dune::KrylovTerminationCriterion::RelativeEnergyError<double> terminationCriterion;
    Dune::Mock::Step step;
  };
}


TEST(RelativeEnergyErrorTerminationCriterion,DeathTests)
{
  auto terminationCriterion = Dune::KrylovTerminationCriterion::RelativeEnergyError<double>();

  ASSERT_DEATH( static_cast<bool>(terminationCriterion), "" );
}

TEST_F(TestRelativeEnergyErrorCriterion, ErrorEstimate)
{
  const auto lookAhead = 5u;
  terminationCriterion.setLookAhead( lookAhead );
  terminationCriterion.setRelativeAccuracy( 1e-3 );

  ASSERT_EQ( terminationCriterion.errorEstimate(), sqrt( std::numeric_limits<double>::max() ) );
  ASSERT_FALSE( static_cast<bool>(terminationCriterion) );

  auto norm = terminationCriterion.relativeAccuracy() / ( lookAhead + 1 );
  step.preconditionedResidualNorm_ = norm*norm;

  for( auto i = 1u; i < lookAhead; ++i )
  {
    ASSERT_EQ( terminationCriterion.errorEstimate(), sqrt( std::numeric_limits<double>::max() ) );
    ASSERT_FALSE( static_cast<bool>(terminationCriterion) );
  }

  ASSERT_TRUE( static_cast<bool>(terminationCriterion) );

  auto denom = lookAhead * norm * norm;
  auto div = 1 + denom;
  ASSERT_DOUBLE_EQ( terminationCriterion.errorEstimate(), sqrt( denom / div ) );
}

