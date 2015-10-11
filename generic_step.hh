#ifndef DUNE_GENERIC_STEP_HH
#define DUNE_GENERIC_STEP_HH

#include <string>
#include <utility>

#include "Mixins/eps.hh"
#include "Mixins/verbosity.hh"
#include "util.hh"

namespace Dune
{
  /*! @cond */
  namespace GenericStepDetail
  {
    class NoData
    {};

    template <class>
    class NoInterface
    {};

    class Ignore
    {
    public:
      template <class... Args>
      void operator()(Args&&...){}
    };
  }
  /*! @endcond */


  /*!
    @ingroup ISTL_Solvers
    @brief Generic policy-based step of an iterative method.

    Solves a linear operator equation \f$Ax=b\f$ with \f$A:\ X\mapsto Y\f$, resp. one of its preconditioned versions
    \f$PAx=Px\f$ or \f$P_1AP_2P_2^{-1}x=P_1b\f$.

      The following steps are performed and may be adjusted:
       1. Apply preconditioner
       2. Compute search direction
       3. Compute scaling for the search direction
       4. Update iterate

    @tparam Domain type of the domain space \f$X\f$
    @tparam Range type of the range space \f$Y\f$
   */
  template <class Domain, class Range,
            class ApplyPreconditioner = GenericStepDetail::Ignore,
            class SearchDirection = GenericStepDetail::Ignore ,
            class Scaling = GenericStepDetail::Ignore ,
            class TreatNonconvexity = GenericStepDetail::Ignore,
            class UpdateIterate = GenericStepDetail::Ignore ,
            class AdjustData = GenericStepDetail::Ignore,
            class Data = GenericStepDetail::NoData ,
            template <class> class Interface = GenericStepDetail::NoInterface>
  class GenericStep :
      public Interface<Data>,
      public Mixin::Eps< real_t<Domain> >,
      public Mixin::Verbosity
  {
    using Interface<Data>::data_;
  public:
    //! type of the domain space
    using domain_type = Domain;
    //! type of the range space
    using range_type = Range;
    //! underlying field type
    using field_type = typename Domain::field_type;
    //! corresponding real type (same as real type for real spaces, differs for complex spaces)
    using real_type = typename FieldTraits<field_type>::real_type;

    template <class... Args>
    GenericStep(Args&&... args)
      : Interface<Data>(std::forward<Args>(args)...)
    {}

    /*!
      @brief Initialize conjugate gradient step.

      This method must be called before starting the conjugate gradient iteration.
      It is responsible for initializing the relevant quantities.

      @param x initial iterate
      @param b initial right hand side
     */
    void init(domain_type& x, range_type& b)
    {
      applyPreconditioner_.pre(data_,x,b);
      data_.init(x,b);
    }

    void reset(domain_type& x, range_type& b)
    {
      data_.reset(x,b);
    }

    /*!
      @brief Perform one step of an iterative method.

      @param x current iterate
      @param b current right hand side
     */
    void compute(domain_type& x, range_type& b)
    {
      applyPreconditioner_(data_);
      computeSearchDirection_(data_);
      computeScaling_(data_);
      treatNonconvexity_(data_,x,this->verbosityLevel());
      updateIterate_(data_,x);
      adjustData_(data_,x);
    }

    /*!
      @brief Post-process final iterate, i.e. apply @code{.cpp} P_.post(x) @endcode
      @param x final iterate
     */
    void postProcess(domain_type& x)
    {
      applyPreconditioner_.post(data_,x);
    }

  private:
    ApplyPreconditioner applyPreconditioner_;
    SearchDirection computeSearchDirection_;
    Scaling computeScaling_;
    TreatNonconvexity treatNonconvexity_;
    UpdateIterate updateIterate_;
    AdjustData adjustData_;
  };
}

#endif // DUNE_GENERIC_STEP_HH
