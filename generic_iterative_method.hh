#ifndef DUNE_GENERIC_ITERATIVE_METHOD_HH
#define DUNE_GENERIC_ITERATIVE_METHOD_HH

#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include <dune/common/typetraits.hh>

#include "optional.hh"
#include "mixins.hh"

#include "fglue/TMP/bind.hh"
#include "fglue/TMP/createMissingBaseClasses.hh"

namespace Dune
{
  //! @cond
  namespace Detail
  {
    using namespace FGlue;

    template <class Step>
    using EnableVerbosity = Apply< StoreIf< IsNotDerivedFrom<Step> > , Mixin::Verbosity >;

    template <class Step,
              class TerminationCriterion,
              class real_type = real_t<typename Step::domain_type> >
    using EnableAdditionalMixinsFromTerminationCriterion =
    Apply< Variadic< StoreIf<
      Apply< Delay<And>,
        FGlue::IsBaseOf<TerminationCriterion>,
        IsNotBaseOf<Step>
        > > ,
      Compose>,
      DUNE_ISTL_MIXINS(real_type)
    >;

    template <class Step,
              class TerminationCriterion>
    using AddMixins =
    Apply< Compose,
      EnableAdditionalMixinsFromTerminationCriterion<Step,TerminationCriterion>,
      EnableVerbosity<Step>
    >;
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
      public Step_ ,
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
      : Step(std::move(step)) ,
        Mixin::MaxSteps(maxSteps) ,
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

    GenericIterativeMethod(const GenericIterativeMethod& other)
      : Step(static_cast<const Step&>(other)),
        Mixin::MaxSteps(other),
        terminate_(other.terminate_)
    {
      initializeConnections();
    }

    GenericIterativeMethod(GenericIterativeMethod&& other)
      : Step(std::move(static_cast<Step&&>(other))),
        Mixin::MaxSteps(std::move(other)),
        terminate_(std::move(other.terminate_))
    {
      initializeConnections();
    }

    GenericIterativeMethod& operator=(const GenericIterativeMethod& other)
    {
      Step::operator=(static_cast<const Step&>(other));
      Mixin::MaxSteps::operator=(other);
      terminate_ = other.terminate_;
      initializeConnections();
    }

    GenericIterativeMethod& operator=(GenericIterativeMethod&& other)
    {
      Step::operator=(std::move(static_cast<Step&&>(other)));
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
      if( this->verbosityLevel() > 1) std::cout << "\n === " << Step::name() << " === " << std::endl;

      initialize(x,b);

      auto step=1u;
      real_type lastErrorEstimate = 1;

      for(; step<=maxSteps(); ++step)
      {
        Step::compute(x,b);

        if( terminate_ )
          break;

        if( Optional::restart(*this) )
        {
          restoreInitialInput(x,b);
          Step::reset(x,b);
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

      Step::postProcess(x);
//      terminate_.print(res);
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
    virtual void apply(domain_type &x, range_type &b, double relativeAccuracy, InverseOperatorResult &res)
    {
      terminate_.setRelativeAccuracy(relativeAccuracy);
      apply(x,b,res);
    }

    /*!
      @brief Apply iterative method to solve \f$Ax=b\f$.
      @param x initial iterate
      @param b initial right hand side
     */
    void apply(domain_type &x, range_type &b)
    {
      InverseOperatorResult res;
      apply(x,b,res);
    }

    //! Access termination criterion.
    TerminationCriterion& getTerminationCriterion()
    {
      return terminate_;
    }

  private:
    /// Initialize connections between iterative method and termination criterion.
    void initializeConnections()
    {
      // connect termination criterion to step implementation to access relevant data
      terminate_.connect(*this);
//      Optional::bind_connect_minimalDecreaseAchieved(terminate_,*this);

      // attach mixins to correctly forward parameters to the termination criterion
      using namespace Mixin;
      Optional::Mixin::Attach< DUNE_ISTL_MIXINS(real_type) >::apply(*this,terminate_);
    }

    void initialize(domain_type& x, range_type& b)
    {
      storeInitialInput(x,b);
      Step::init(x,b);
      terminate_.init();
    }

    void storeInitialInput(const domain_type& x, const range_type& b)
    {
      x0 = std::unique_ptr<domain_type>(new domain_type(x));
      b0 = std::unique_ptr<range_type>(new range_type(b));
    }

    void restoreInitialInput(domain_type& x, range_type& b)
    {
      x = *x0;
      b = *b0;
    }

    void printOutput(unsigned step, real_type lastErrorEstimate) const
    {
      this->printHeader(std::cout);
      InverseOperator<domain_type,range_type>::printOutput(std::cout,
                                                           static_cast<decltype(terminate_.errorEstimate())>(step),
                                                           terminate_.errorEstimate(),
                                                           lastErrorEstimate);
    }

    void printFinalOutput(const InverseOperatorResult& res, unsigned step) const
    {
      auto name = this->name();
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

    TerminationCriterion terminate_;
    std::unique_ptr<domain_type> x0;
    std::unique_ptr<range_type> b0;
  };

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
