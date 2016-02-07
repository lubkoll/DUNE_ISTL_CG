#ifndef DUNE_GENERIC_ITERATIVE_METHOD_HH
#define DUNE_GENERIC_ITERATIVE_METHOD_HH

#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include "dune/common/typetraits.hh"
#include "dune/istl/solver.hh"

#include "optional.hh"
#include "mixins.hh"

#include "fglue/TMP/bind.hh"
#include "fglue/TMP/createMissingBaseClasses.hh"

namespace Dune
{
  //! @cond
  namespace Detail
  {
    /// Empty default storage object.
    template <class domain_type, class range_type, class TerminationCriterion, class = void>
    struct Storage
    {
      void store(const domain_type&, const range_type&) const noexcept
      {}

      void restore(domain_type&, range_type&) const noexcept
      {}
    };

    /**
     * @brief Storage object for GenericIterativeMethod with a termination criterion that may trigger restarts.
     *
     * Stores initial guess \f$ x0 \f$ and initial right hand side \f$ b0 \f$.
     */
    template <class domain_type, class range_type, class TerminationCriterion>
    struct Storage< domain_type, range_type, TerminationCriterion, void_t< Try::MemFn_restart<TerminationCriterion> > >
    {
      Storage(Storage&&) = default;
      Storage& operator=(Storage&&) = default;

      Storage(const Storage& other)
        : x0( other.x0 ? new domain_type(*other.x0) : nullptr ),
          b0( other.b0 ? new range_type(*other.b0) : nullptr)
      {}

      Storage& operator=(const Storage& other)
      {
        if( other.x0 )
          x0.reset( new domain_type(*other.x0) );
        if( other.b0 )
          b0.reset( new range_type(*other.b0) );
      }

      void store(const domain_type& x, const range_type& b)
      {
        x0.reset( new domain_type(x) );
        b0.reset( new range_type(b) );
      }

      void restore(domain_type& x, range_type& b) const
      {
        assert(x0 && b0);
        x = *x0;
        b = *b0;
      }

      std::unique_ptr<domain_type> x0;
      std::unique_ptr<range_type> b0;
    };


    using namespace FGlue;

    /// Is Empty if Step is derived from Mixin::Verbosity, else is Mixin::Verbosity.
    template <class Step>
    using EnableVerbosity = Apply< StoreIf< IsNotDerivedFrom<Step> > , Mixin::Verbosity >;

    /**
     * Generate a template meta-function that takes an arbitrary number of arguments and generates a
     * type that is derived from each argument that is a base class of TerminationCriterion, but not
     * of Step.
     */
    template <class Step,
              class TerminationCriterion>
    using AdditionalMixinsCondition =
    Apply< Delay<Or>,
      FGlue::IsBaseOf<TerminationCriterion>,
      FGlue::IsBaseOf<Step>
    >;

    /// Generate type that is derived from all necessary mixin base classes for GenericIterativeMethod.
    template <class Step,
              class TerminationCriterion,
              class real_type = real_t<typename Step::domain_type> >
    using AddMixins =
    Apply< Compose,
      EnableBaseClassesIf< AdditionalMixinsCondition<Step,TerminationCriterion> , DUNE_ISTL_MIXINS( real_type ) >,
      EnableVerbosity<Step>
    >;
  }

  namespace Optional
  {
    template < class Step >
    using TryNestedType_Cache = typename Step::Cache;

    struct NoCache
    {
      template <class... Args>
      explicit NoCache(Args&&...) {}
    };

    template < class Step , class = void >
    struct StepTraits
    {
      using Cache = NoCache;

      static void setCache( const Step&, Cache* ) noexcept
      {}
    };

    template < class Step >
    struct StepTraits< Step, void_t< TryNestedType_Cache<Step> > >
    {
      using Cache = TryNestedType_Cache<Step>;

      static void setCache( Step& step, Cache* cache )
      {
        step.setCache(cache);
      }
    };

    template < class Step, class domain_type, class range_type >
    typename StepTraits< Step >::Cache createCache( domain_type& x, range_type& y )
    {
      return typename StepTraits< Step >::Cache( x, y );
    }

    template < class Step, class Cache >
    void setCache( Step& step, Cache* cache )
    {
      StepTraits< Step >::setCache( step, cache );
    }
  }
  //! @endcond


  /*!
    @ingroup ISTL_Solvers
    @brief Generic wrapper for iterative methods.
   */
  template <class Step_,
            class TerminationCriterion_,
            class real_type = real_t<typename Step_::domain_type> >
  class GenericIterativeMethod :
      public InverseOperator<typename Step_::domain_type, typename Step_::range_type> ,
      public Mixin::MaxSteps ,
      public Detail::AddMixins<Step_,TerminationCriterion_>
  {
  public:
    using Step = Step_;
    using TerminationCriterion = TerminationCriterion_;
    using domain_type = typename Step::domain_type;
    using range_type  = typename Step::range_type;
    using field_type  = field_t<domain_type>;

    /*!
      @brief Construct from given step implementation and termination criterions.
      @param step object implementing one step of an iterative scheme
      @param terminate termination criterion
      @param maxSteps
     */
    GenericIterativeMethod(Step step, TerminationCriterion terminate, unsigned maxSteps = 1000)
      : Mixin::MaxSteps(maxSteps) ,
        step_(std::move(step)) ,
        terminate_(std::move(terminate))
    {
      initializeConnections();
    }

    //! @brief Optional constructor for the case that the step has a constructor that takes three parameters of type Operator, Preconditioner and ScalarProduct.
    template <class Operator, class Preconditioner, class ScalarProduct,
              typename std::enable_if<std::is_constructible<Step,Operator,Preconditioner,ScalarProduct>::value>::type* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, ScalarProduct&& sp, TerminationCriterion terminate, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P),std::forward<ScalarProduct>(sp)) ,
                                std::move(terminate) ,
                                maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes three parameters of type Operator, Preconditioner and ScalarProduct and the termination criterion is default constructible.
    template <class Operator, class Preconditioner, class ScalarProduct,
              typename std::enable_if<std::is_constructible<Step,Operator,Preconditioner,ScalarProduct>::value && std::is_default_constructible<TerminationCriterion>::value>::type* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, ScalarProduct&& sp, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P),std::forward<ScalarProduct>(sp)) ,
                                TerminationCriterion() ,
                                maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes three parameters of type Operator and Preconditioner.
    template <class Operator, class Preconditioner,
              typename std::enable_if<std::is_constructible<Step,Operator,Preconditioner>::value>::type* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, TerminationCriterion terminate, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P)) ,
                                std::move(terminate) ,
                                maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes two parameters of type Operator and Preconditioner and the termination criterion is default constructible.
    template <class Operator, class Preconditioner,
              typename std::enable_if<std::is_constructible<Step,Operator,Preconditioner>::value && std::is_default_constructible<TerminationCriterion>::value>::type* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P)) ,
                                TerminationCriterion() ,
                                maxSteps )
    {}

//    GenericIterativeMethod(const GenericIterativeMethod& other)
//      : Step(static_cast<const Step&>(other)),
//        Mixin::MaxSteps(other),
//        terminate_(other.terminate_)
//    {
//      initializeConnections();
//    }

    GenericIterativeMethod(GenericIterativeMethod&& other)
      : Mixin::MaxSteps( other.maxSteps() ),
        step_( std::move( other.step_ ) ),
        terminate_( std::move( other.terminate_ ) )
    {
      initializeConnections();
    }

//    GenericIterativeMethod& operator=(const GenericIterativeMethod& other)
//    {
//      Step::operator=(static_cast<const Step&>(other));
//      Mixin::MaxSteps::operator=(other);
//      terminate_ = other.terminate_;
//      initializeConnections();
//    }

    GenericIterativeMethod& operator=(GenericIterativeMethod&& other)
    {
      step_ = std::move(static_cast<Step&&>(other));
      Mixin::MaxSteps::operator=(std::move(other));
      terminate_ = std::move(other.terminate_);
      initializeConnections();
    }

    /*!
      @brief Apply iterative method to solve \f$Ax=b\f$.
      @param x initial iterate
      @param b initial right hand side
      @param res some statistics
     */
    virtual void apply(domain_type& x, range_type& b, InverseOperatorResult& res)
    {
      if( this->verbosityLevel() > 1)
        std::cout << "\n === " << step_.name() << " === " << std::endl;

      auto cache = Optional::createCache< Step >( x, b );
      Optional::setCache( step_, &cache );

      initialize(x,b);

      auto step=1u;
      real_type lastErrorEstimate = 1;

      for(; step<=maxSteps(); ++step)
      {
        step_.compute(x,b);

        if( terminate_ )
          break;

        if( Optional::restart( step_ ) )
        {
          storage_.restore(x,b);
          step_.reset(x,b);
          terminate_.init();
          step = 0u;
          lastErrorEstimate = 1;
          continue;
        }

        if( this->verbosityLevel() > 1 )
        {
          printOutput(step,lastErrorEstimate);
          lastErrorEstimate = terminate_.errorEstimate();
        }
      }

      step_.postProcess(x);
      terminate_.print(res);
      if( step < maxSteps() + 1 ) res.converged = true;
      if( this->is_verbose() )  printFinalOutput(res,step);
    }

    /*!
      @brief Apply iterative method to solve \f$Ax=b\f$.
      @param x initial iterate
      @param b initial right hand side
      @param relativeAccuracy required relative accuracy
      @param res some statistics
     */
    virtual void apply( domain_type &x, range_type &b, double relativeAccuracy, InverseOperatorResult &res )
    {
      terminate_.setRelativeAccuracy( relativeAccuracy );
      apply( x, b, res );
    }

    /*!
      @brief Apply iterative method to solve \f$Ax=b\f$.
      @param x initial iterate
      @param b initial right hand side
     */
    void apply( domain_type &x, range_type &b )
    {
      InverseOperatorResult res;
      apply( x, b, res);
    }

    //! Access termination criterion.
    TerminationCriterion& getTerminationCriterion()
    {
      return terminate_;
    }

    //! Access step implementation.
    Step& getStep()
    {
      return step_;
    }

  private:
    /// Initialize connections between iterative method and termination criterion.
    void initializeConnections()
    {
      // connect termination criterion to step implementation to access relevant data
      terminate_.connect( step_ );
//      Optional::bind_connect_minimalDecreaseAchieved(terminate_,*this);

      // attach mixins to correctly forward parameters to the termination criterion
      using namespace Mixin;
      Optional::Mixin::Attach< DUNE_ISTL_MIXINS( real_type ) >::apply( *this, step_ );
      Optional::Mixin::Attach< DUNE_ISTL_MIXINS( real_type ) >::apply( *this, terminate_ );
    }

    void initialize(domain_type& x, range_type& b)
    {
      storage_.store(x,b);
      step_.init(x,b);
      terminate_.init();
    }

    void printOutput( unsigned step, real_type lastErrorEstimate ) const
    {
      this->printHeader( std::cout );
      InverseOperator< domain_type, range_type >::printOutput( std::cout,
                                                               static_cast<decltype(terminate_.errorEstimate())>(step),
                                                               terminate_.errorEstimate(),
                                                               lastErrorEstimate );
    }

    void printFinalOutput(const InverseOperatorResult& res, unsigned step) const
    {
      auto name = step_.name();
      name += ( (step==maxSteps()+1) ? ": Failed" : ": Converged" );
      std::cout << "\n === " << name << " === " << std::endl;
      this->printHeader(std::cout);
      InverseOperator<domain_type,range_type>::printOutput(std::cout,
                                                           static_cast<decltype(res.reduction)>(res.iterations),
                                                           res.reduction,
                                                           res.reduction/res.conv_rate);
      name = std::string(name.size(),'=');
      std::cout << " === " << name << "" << " === \n" << std::endl;
    }

    Step step_;
    TerminationCriterion terminate_;
    Detail::Storage<domain_type,range_type,TerminationCriterion> storage_;
  };

  /*!
    @brief Generating function for GenericIterativeMethod.

    @param step Implementation of one step of an iterative method.
    @param terminationCriterion Implementation of a termination criterion.
   */
  template <class Step, class TerminationCriterion>
  GenericIterativeMethod< typename std::decay<Step>::type, typename std::decay<TerminationCriterion>::type >
  makeGenericIterativeMethod(Step&& step, TerminationCriterion&& terminationCriterion)
  {
    return GenericIterativeMethod< typename std::decay<Step>::type, typename std::decay<TerminationCriterion>::type >
        ( std::forward<Step>(step),
          std::forward<TerminationCriterion>(terminationCriterion) );
  }
}

#endif // DUNE_GENERIC_ITERATIVE_METHOD_HH
