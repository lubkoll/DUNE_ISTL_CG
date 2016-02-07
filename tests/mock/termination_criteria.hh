#ifndef DUNE_ISTL_TESTS_MOCK_TERMINATION_CRITERIA_HH
#define DUNE_ISTL_TESTS_MOCK_TERMINATION_CRITERIA_HH

#include "mixins/absoluteAccuracy.hh"
#include "mixins/minimalAccuracy.hh"
#include "mixins/relativeAccuracy.hh"
#include "mixins/eps.hh"
#include "mixins/verbosity.hh"

#include "optional.hh"

namespace Dune
{
  namespace Mock
  {
    template <class Step>
    struct TerminationCriterion : Dune::Mixin::RelativeAccuracy<double>
    {
      TerminationCriterion(bool val = true)
        : value(val)
      {}

      void init()
      {
        wasInitialized = true;
      }

      operator bool() const
      {
        if( step_ != nullptr && Optional::terminate(*step_))
          return true;

        return value;
      }

      void connect(const Step& step)
      {
        step_ = &step;
      }

      double absoluteError() const
      {
        return 1;
      }

      double errorEstimate() const
      {
        return 1;
      }

      template <class InverseOperatorResult>
      void print(InverseOperatorResult& res)
      {}

      bool wasInitialized = false;
      bool value = true;
      const Step* step_ = nullptr;
    };

    template <class Step>
    struct MixinTerminationCriterion
        : Dune::Mixin::AbsoluteAccuracy<double>,
          Dune::Mixin::RelativeAccuracy<double>,
          Dune::Mixin::MinimalAccuracy<double>,
          Dune::Mixin::Eps<double>,
          Dune::Mixin::Verbosity
    {
      void init()
      {
        wasInitialized = true;
      }

      operator bool() const
      {
        return true;
      }

      void connect(const Step&)
      {}

      double absoluteError() const
      {
        return 1;
      }

      double errorEstimate() const
      {
        return 1;
      }

      template <class InverseOperatorResult>
      void print(InverseOperatorResult& res)
      {}

      bool wasInitialized = false;
    };
  }
}

#endif // DUNE_ISTL_TESTS_MOCK_TERMINATION_CRITERIA_HH
