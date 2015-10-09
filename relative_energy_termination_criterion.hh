#ifndef DUNE_TERMINATION_CRITERIA_HH
#define DUNE_TERMINATION_CRITERIA_HH

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include <dune/common/timer.hh>

#include "mixins.hh"
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
      @brief %Termination criterion for conjugate gradient methods based on an estimate of the relative energy error.

      Relative energy error termination criterion according to @cite Strakos2005 (see also @cite Hestenes1952, and for a related absolute energy error
      criterion @cite Arioli2004).

      Requires that CG starts at \f$ x = 0 \f$. More general starting values might be used, but must be chosen such that
      the estimate for the energy norm of the solution stays positive (see the above mentioned paper for details).

      The essential idea behind this termination criterion is simple: perform \f$d\f$ extra iterations of the conjugate gradient
      method to estimate the absolute or relative error in the energy norm (the parameter \f$d\f$ can be adjusted with setLookAhead()).
      To compute the error estimate only quantities that are anyway computed as intermediate results in the conjugate gradient method are required.

      This estimate only relies on local orthogonality and thus its evaluation is numerically stable.
     */
    template <class real_type>
    class RelativeEnergyError :
        public Mixin::AbsoluteAccuracy<real_type>,
        public Mixin::RelativeAccuracy<real_type>,
        public Mixin::MinimalAccuracy<real_type>,
        public Mixin::Eps<real_type>
    {
      /*! @cond */
      class TypeErasedCGHolder
      {
        class AbstractBase
        {
        public:
          virtual ~AbstractBase(){}
          virtual real_type alpha() const = 0;
          virtual real_type length() const = 0;
          virtual real_type preconditionedResidualNorm() const = 0;
        };

        template <class Type>
        class Base : public AbstractBase
        {
        public:
          Base(const Type& type) : type_(&type) {}
          real_type alpha() const final override { return type_->alpha(); }
          real_type length() const final override { return type_->length(); }
          real_type preconditionedResidualNorm() const final override { return type_->preconditionedResidualNorm(); }
        private:
          const Type* type_;
        };

      public:
        TypeErasedCGHolder() = default;

        template <class Type>
        TypeErasedCGHolder& operator=(const Type& type)
        {
          impl_ = std::make_unique< Base<Type> >(type);
          return *this;
        }

        template <class Type,
                  access_t< std::enable_if< std::is_rvalue_reference<Type>::value > >* = nullptr >
        TypeErasedCGHolder& operator=(Type&&)
        {
          throw std::invalid_argument("TypeErasedCGHolder can not be assigned with a temporary object.");
          return *this;
        }

        real_type alpha() const { return impl_->alpha(); }

        real_type length() const { return impl_->length(); }

        real_type preconditionedResidualNorm() const { return impl_->preconditionedResidualNorm(); }

        operator bool() const { return impl_!=nullptr; }

      private:
        std::unique_ptr<AbstractBase> impl_ = nullptr;
      };
      /*! @endcond */
    public:
      /*!
        @brief Constructor.
        @param relativeAccuracy required relative accuracy for the estimated energy error
        @param verbosity verbosity level
        @param eps maximal attainable accuracy
        @param absoluteAccuracy absolute accuracy
       */
      RelativeEnergyError(real_type relativeAccuracy = std::numeric_limits<real_type>::epsilon(),
                          unsigned verbosity = 0,
                          real_type eps = std::numeric_limits<real_type>::epsilon(),
                          real_type absoluteAccuracy = std::numeric_limits<real_type>::epsilon())
        : Mixin::AbsoluteAccuracy<real_type>{absoluteAccuracy},
          Mixin::RelativeAccuracy<real_type>{relativeAccuracy},
          Mixin::Eps<real_type>{eps}
      {}

      //! @copydoc ResidualBased::operator bool()
      operator bool()
      {
        readParameter();

        if( vanishingStep() ) return true;

        using std::max;
        auto acc = max( this->relativeAccuracy() , this->eps() );
        return scaledGamma2.size() > lookAhead_ && errorEstimate() < acc;
      }

      //! @copydoc ResidualBased::init()
      void init()
      {
        scaledGamma2.clear();
        energyNorm2 = stepLength2 = 0;
        watch.reset();
        watch.start();
      }

      /*!
        @brief Access estimated error.
       */
      auto errorEstimate() const
      {
        return sqrt(squaredRelativeError());
      }

      /*!
        @brief Set the additional CG-iterations required for estimating the relative energy error.

        @param lookAhead the requested lookahead value (default = 5)
       */
      void setLookAhead(unsigned lookAhead = 5)
      {
        lookAhead_ = lookAhead;
      }

      /*!
        @brief Relaxed termination criterion.
        @return true if the iteration has reached some minimal required accuracy, possibly bigger than the desired accuracy. This method is required in the
         truncated regularized conjugate gradient method only.
       */
      bool minimalDecreaseAchieved() const
      {
        return squaredRelativeError() < this->minimalAccuracy()*this->minimalAccuracy();
      }


      //! @copydoc ResidualBased::connect()
      template <class Step>
      void connect(Step&& step)
      {
        step_ = std::forward<Step>(step);
      }

      //! @copydoc ResidualBased::print()
      void print(InverseOperatorResult& res)
      {
        res.iterations = scaledGamma2.size();
        res.reduction = sqrt(squaredRelativeError());
        res.conv_rate = pow(res.reduction,1./res.iterations);
        res.elapsed = watch.stop();
      }

      /*!
        @brief check if the energy norm of the current step \f$\|q\|_A=\sqrt(qAq)\f$ is smaller than the maximal attainable accuracy multiplied with the energy norm of the iterate \f$\varepsilon_{max}\|x\|_A\f$.
        @return true if \f$\|q\|<\varepsilon_{max}\|x\|_A\f$, else false
       */
      bool vanishingStep() const
      {
        auto acc2 = this->absoluteAccuracy()*this->absoluteAccuracy();
        using std::min;
        if( energyNorm2 > acc2) acc2 = min(acc2,this->eps()*this->eps()*energyNorm2);
        return stepLength2 < acc2;
      }

    private:
      void readParameter()
      {
        assert(step_);
        scaledGamma2.push_back( step_.alpha() * step_.preconditionedResidualNorm() );
        energyNorm2 += scaledGamma2.back();
        using std::abs;
        stepLength2 = abs(step_.length());
      }

      real_type squaredRelativeError() const
      {
        if( scaledGamma2.size() < lookAhead_ ) return std::numeric_limits<real_type>::max();
        return std::accumulate(scaledGamma2.end() - lookAhead_, scaledGamma2.end(), real_type(0)) / energyNorm2;
      }

      unsigned lookAhead_ = 25;
      std::vector<real_type> scaledGamma2 = std::vector<real_type>{};
      real_type energyNorm2 = 0;
      real_type stepLength2 = 0;
      TypeErasedCGHolder step_ = {};
      Timer watch = Timer{false};
    };
  }
}

#endif // DUNE_TERMINATION_CRITERIA_HH
