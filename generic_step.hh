#ifndef DUNE_GENERIC_STEP_HH
#define DUNE_GENERIC_STEP_HH

#include <string>
#include <utility>

#include "mixins.hh"
#include "tmp/compile_time_sequence.hh"
#include "tmp/logic.hh"
#include "tmp/optional_base_class.hh"
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
    @brief Generic step of an iterative method.

    Solves a linear operator equation \f$Ax=b\f$ with \f$A:\ X\mapsto Y\f$, resp. one of its preconditioned versions
    \f$PAx=Px\f$ or \f$P_1AP_2P_2^{-1}x=P_1b\f$.

      The following steps are performed and may be adjusted:
       1. Apply preconditioner
       2. Compute search direction
       3. Compute scaling for the search direction
       4. Possibly treat nonconvexity
       5. Update iterate
       6. Adjust other internal data (such as the residual)

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
      public TMP::BaseClassesIf<
        TMP::OrUnaryToSequence<
          TMP::IsDerivedFrom ,
          TMP::Sequence<ApplyPreconditioner,SearchDirection,Scaling,TreatNonconvexity,UpdateIterate,AdjustData,Data>
        >,
        Mixin::IterativeRefinements , Mixin::Verbosity , Mixin::Eps< real_t<Domain> >
      >::type
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
    {
      initializeConnections();
    }

    GenericStep(const GenericStep& other)
      : Interface<Data>(other.data_)
    {
      initializeConnections();
    }

    GenericStep& operator=(const GenericStep& other)
    {
      Interface<Data>::operator=(other.data_);
      initializeConnections();
    }

    GenericStep& operator=(GenericStep&& other)
    {
      Interface<Data>::operator=(std::move(other.data_));
      initializeConnections();
    }

    GenericStep(GenericStep&& other)
      : Interface<Data>(std::move(other.data_))
    {
      initializeConnections();
    }

    /*!
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
      treatNonconvexity_(data_,x);
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
    void initializeConnections()
    {
      using namespace Mixin;
      Optional::Mixin::Attach< IterativeRefinements , Verbosity , Eps<real_type> >::apply(*this,
                                                                                          applyPreconditioner_,
                                                                                          computeSearchDirection_,
                                                                                          computeScaling_,
                                                                                          treatNonconvexity_,
                                                                                          updateIterate_,
                                                                                          adjustData_,
                                                                                          data_);
    }

    ApplyPreconditioner applyPreconditioner_;
    SearchDirection computeSearchDirection_;
    Scaling computeScaling_;
    TreatNonconvexity treatNonconvexity_;
    UpdateIterate updateIterate_;
    AdjustData adjustData_;
  };
}

#endif // DUNE_GENERIC_STEP_HH
