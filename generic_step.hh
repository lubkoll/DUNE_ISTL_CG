#ifndef DUNE_GENERIC_STEP_HH
#define DUNE_GENERIC_STEP_HH

#include <string>
#include <utility>

#include <dune/common/typetraits.hh>

#include "mixins.hh"
#include "fglue/TMP/createMissingBaseClasses.hh"
#include "fglue/Fusion/connect.hh"

namespace Dune
{
  /*! @cond */
  namespace GenericStepDetail
  {
    class NoCache
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


    template < class ApplyPreconditioner,
               class ComputeSearchDirection,
               class ComputeScaling,
               class Update,
               class real_type >
    using AddMixins =
    FGlue::EnableBaseClassesIf<
      FGlue::IsBaseOfOneOf< ApplyPreconditioner, ComputeSearchDirection, ComputeScaling, Update >,
      DUNE_ISTL_MIXINS( real_type )
    >;

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
            class ApplyPreconditioner    = GenericStepDetail::Ignore,
            class ComputeSearchDirection = GenericStepDetail::Ignore,
            class ComputeScaling         = GenericStepDetail::Ignore,
            class Update                 = GenericStepDetail::Ignore,
            class Interface              = GenericStepDetail::Ignore>
  class GenericStep :
      public GenericStepDetail::AddMixins< ApplyPreconditioner, ComputeSearchDirection, ComputeScaling, Update, real_t<Domain> >,
      public Interface
  {
  public:
    //! type of the domain space
    using domain_type = Domain;
    //! type of the range space
    using range_type = Range;
    //! underlying field type
    using field_type = field_t<Domain>;
    //! corresponding real type (same as real type for real spaces, differs for complex spaces)
    using real_type = real_t<Domain>;
    //! cache object storing temporaries
    using Cache = typename Interface::Cache;

    template <class LinOp, class Prec, class SP>
    GenericStep(LinOp& A, Prec& P, SP& sp)
      : A_(A), P_(P), ssp_(), sp_(sp)
    {
      initializeConnections();
//        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
//        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
    }

    template <class LinOp, class Prec>
    GenericStep(LinOp& A,  Prec& P)
      : A_(A), P_(P), ssp_(), sp_(ssp_)
    {
      initializeConnections();
//        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
//        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
    }

    GenericStep( const GenericStep& other )
      : A_( other.A_ ),
        P_( other.P_ ),
        ssp_( ),
        sp_( other.sp_ )
    {
      initializeConnections( );
    }

    GenericStep( GenericStep&& other )
      : A_( other.A_ ),
        P_( other.P_ ),
        ssp_(),
        sp_( other.sp_ )
    {
      initializeConnections( );
    }

    /*!
      @param x initial iterate
      @param b initial right hand side
     */
    void init(domain_type& x, range_type& b)
    {
      applyPreconditioner_.pre( P_, x, b );
    }

    void reset(domain_type& x, range_type& b)
    {
      this->cache_->reset( &A_, &P_, &sp_ );
    }

    void setCache( Cache* cache )
    {
      Interface::setCache( cache );
      this->cache_->reset( &A_, &P_, &sp_ );
    }

    /*!
      @brief Perform one step of an iterative method.

      @param x current iterate
      @param b current right hand side
     */
    void compute(domain_type&, range_type&)
    {
      applyPreconditioner_( *this->cache_ );
      computeSearchDirection_( *this->cache_ );
      computeScaling_( *this->cache_ );
      update_( *this->cache_ );
    }



    /*!
      @brief Post-process final iterate, i.e. apply @code{.cpp} P_.post(x) @endcode
      @param x final iterate
     */
    void postProcess(domain_type& x)
    {
      applyPreconditioner_.post(P_,x);
    }

  private:
    void initializeConnections()
    {
      using FGlue::Connector;
      using FGlue::IsDerivedFrom;

      Connector< IsDerivedFrom<Mixin::IterativeRefinements> >::template
          from< Mixin::IterativeRefinements >( *this ).
          to(applyPreconditioner_,computeSearchDirection_,computeScaling_,update_);
      Connector< IsDerivedFrom<Mixin::Verbosity> >::template
          from< Mixin::Verbosity >( *this ).
          to(applyPreconditioner_,computeSearchDirection_,computeScaling_,update_);
      Connector< IsDerivedFrom< Mixin::Eps<real_type> > >::template
          from< Mixin::Eps<real_type> >( *this ).
          to(applyPreconditioner_,computeSearchDirection_,computeScaling_,update_);
    }


    LinearOperator<Domain,Range>& A_;
    Preconditioner<Domain,Range>& P_;
    SeqScalarProduct<Domain> ssp_;
    ScalarProduct<Domain>& sp_;

    ApplyPreconditioner applyPreconditioner_;
    ComputeSearchDirection computeSearchDirection_;
    ComputeScaling computeScaling_;
    Update update_;
  };
}

#endif // DUNE_GENERIC_STEP_HH
