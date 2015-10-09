#ifndef DUNE_RESIDUAL_BASED_TERMINATION_CRITERION_HH
#define DUNE_RESIDUAL_BASED_TERMINATION_CRITERION_HH

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include <dune/common/timer.hh>

#include "Mixins/eps.hh"
#include "Mixins/relativeAccuracy.hh"
#include "util.hh"

namespace Dune
{
  /*! @cond */
  class InverseOperatorResult;
  /*! @endcond */

  namespace KrylovTerminationCriterion
  {
    /*!
      @ingroup ISTL_Solvers
      @brief Residual-based relative error criterion.
     */
    template <class real_type>
    class ResidualBased :
        public Mixin::Eps<real_type> ,
        public Mixin::RelativeAccuracy<real_type>
    {
    public:
      /*!
        @brief Constructor.
        @param accuracy required relative accuracy of the residual
        @param eps maximal attainable accuracy
       */
      ResidualBased(real_type accuracy = std::numeric_limits<real_type>::epsilon(),
                    real_type eps = std::numeric_limits<real_type>::epsilon() )
        : Mixin::Eps<real_type>{eps},
          Mixin::RelativeAccuracy<real_type>{accuracy}
      {}

      /*!
        @brief Initialize internal state before using the termination criterion.
       */
      void init()
      {
        assert(step_residualNorm_);
        initialResidualNorm_ = step_residualNorm_();
        iteration_ = 0;
        watch.reset();
        watch.start();
      }

      /*!
        @brief Connect to step implementation.

        @warning operator bool() seg-faults if no step is connected

        @param step step implementation that provides a member function step.residualNorm()
       */
      template <class Step,
                class = std::enable_if<std::is_reference<Step>::value> >
      void connect(Step&& step)
      {
        step_residualNorm_ = std::bind(&access_t< std::decay<Step> >::residualNorm,&step);
      }

      /*!
        @brief Write information to res.
        @param res holds information on required iterations, reduction, convergence, rate and elapsed time.
       */
      void print(InverseOperatorResult& res)
      {
        assert(step_residualNorm_);
        res.iterations = iteration_;
        res.reduction = step_residualNorm_()/initialResidualNorm_;
        res.conv_rate = pow(res.reduction,1./res.iterations);
        res.elapsed = watch.stop();
      }

      /*!
        @brief Evaluate termination criterion.
        @return true if termination criterion is satisfied, else false
       */
      operator bool()
      {
        ++iteration_;

        auto acc = std::max(this->eps(),this->relativeAccuracy());

        if( errorEstimate() < acc )
          return true;

        return false;
      }

      /*!
        @brief Access estimated error.
       */
      auto errorEstimate() const
      {
        assert(step_residualNorm_);
        return step_residualNorm_()/initialResidualNorm_;
      }

    private:
      real_type initialResidualNorm_ = -1;
      unsigned iteration_ = 0;
      std::function<real_type()> step_residualNorm_;
      Timer watch = Timer{false};
    };
  }
}

#endif // DUNE_RESIDUAL_BASED_TERMINATION_CRITERION_HH
