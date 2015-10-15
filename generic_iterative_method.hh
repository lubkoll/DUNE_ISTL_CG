#ifndef DUNE_GENERIC_ITERATIVE_METHOD_HH
#define DUNE_GENERIC_ITERATIVE_METHOD_HH

#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include "optional.hh"
#include "tmp/for_each.hh"
#include "tmp/logic.hh"
#include "voider.hh"
#include "util.hh"
#include "mixins.hh"

namespace Dune
{
  //! @cond
  class InverseOperatorResult;
  //! @endcond

  inline std::ostream& operator<<(std::ostream& os, const InverseOperatorResult& res)
  {
    os << "InverseOperatorResult: converged         : " << std::boolalpha << res.converged << "\n";
    os << "                       iterations        : " << res.iterations << "\n";
    os << "                       reduction         : " << res.reduction << "\n";
    os << "                       convergence rate  : " << res.conv_rate << "\n";
    os << "                       elapsed time      : " << res.elapsed << " seconds\n";
    return os;
  }

  //! @cond
  namespace GenericIterativeMethodDetail
  {
    using namespace TMP;

    template <class Step>
    using EnableVerbosity = Bind< StoreIfNotDerivedFrom<Step> , Mixin::Verbosity >;

    template <class Step, class TerminationCriterion>
    using MixinOperation = StoreIf< typename Apply< And , NotBaseOf<Step> , BaseOf<TerminationCriterion> >::type >;

    template <class Step, class TerminationCriterion, class... Mixins>
    using EnableMixins = Bind< VariadicApply< MixinOperation<Step,TerminationCriterion> , Compose > , Mixins... >;

    template <class Step, class TerminationCriterion, class... Mixins>
    using Base = typename Apply<
      Compose,
      EnableVerbosity<Step>,
      EnableMixins<Step,TerminationCriterion,Mixins...>
    >::type::template apply<>::type;
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
      public GenericIterativeMethodDetail::Base<
        Step_,
        TerminationCriterion_,
        Mixin::RelativeAccuracy<real_type>, Mixin::AbsoluteAccuracy<real_type>, Mixin::MinimalAccuracy<real_type>, Mixin::Eps<real_type>
      >
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
              std::enable_if_t<std::is_constructible<Step,Operator,Preconditioner,ScalarProduct>::value>* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, ScalarProduct&& sp, TerminationCriterion terminate, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P),std::forward<ScalarProduct>(sp)) ,
                                std::move(terminate) , maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes three parameters of type Operator, Preconditioner and ScalarProduct and the termination criterion is default constructible.
    template <class Operator, class Preconditioner, class ScalarProduct,
              std::enable_if_t<std::is_constructible<Step,Operator,Preconditioner,ScalarProduct>::value && std::is_default_constructible<TerminationCriterion>::value>* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, ScalarProduct&& sp, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P),std::forward<ScalarProduct>(sp)) ,
                                TerminationCriterion() , maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes three parameters of type Operator and Preconditioner.
    template <class Operator, class Preconditioner,
              std::enable_if_t<std::is_constructible<Step,Operator,Preconditioner>::value>* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, TerminationCriterion terminate, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P)) ,
                                std::move(terminate) , maxSteps )
    {}

    //! @brief Optional constructor for the case that the step has a constructor that takes two parameters of type Operator and Preconditioner and the termination criterion is default constructible.
    template <class Operator, class Preconditioner,
              std::enable_if_t<std::is_constructible<Step,Operator,Preconditioner>::value && std::is_default_constructible<TerminationCriterion>::value>* = nullptr >
    GenericIterativeMethod(Operator&& A, Preconditioner&& P, unsigned maxSteps = 1000)
      : GenericIterativeMethod( Step(std::forward<Operator>(A),std::forward<Preconditioner>(P)) ,
                                TerminationCriterion() , maxSteps )
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

      storeInitialInput(x,b);
      Step::init(x,b);
      terminate_.init();

      auto step=1u;
      real_type lastErrorEstimate = 1;
      for(; step<=maxSteps(); ++step)
      {
        computeStep(x,b);

        if( terminate_ || Optional::terminate(*this) ) break;
        if( Optional::restart(*this) )
        {
          restoreInitialInput(x,b);
          Step::reset(x,b);
          terminate_.init();
          step = 1u;
          lastErrorEstimate = 1;
          continue;
        }

        if( this->verbosityLevel() > 1 ) printOutput(step,lastErrorEstimate);
        lastErrorEstimate = terminate_.errorEstimate();
      }

      Step::postProcess(x);
      terminate_.print(res);
      if( step < maxSteps() + 1 ) res.converged = true;
      if( this->verbose() )  printFinalOutput(res,step);
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
    TerminationCriterion& terminationCriterion()
    {
      return terminate_;
    }

  private:
    void initializeConnections()
    {
      terminate_.connect(*this);
      Optional::bind_connect_minimalDecreaseAchieved(terminate_,*this);

      using namespace Mixin;
      Optional::Mixin::Attach< AbsoluteAccuracy<real_type>, MinimalAccuracy<real_type>, RelativeAccuracy<real_type> , Verbosity , Eps<real_type> >::apply(*this,terminate_);
    }

    void storeInitialInput(const domain_type& x, const range_type& b)
    {
      x0 = std::make_unique<domain_type>(x);
      b0 = std::make_unique<range_type>(b);
    }

    void restoreInitialInput(domain_type& x, range_type& b)
    {
      x = *x0;
      b = *b0;
    }

    void computeStep(domain_type& x, range_type& b)
    {
      try { Step::compute(x,b); }
      catch(...)
      {
        restoreInitialInput(x,b);
        throw;
      }
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
}

#endif // DUNE_GENERIC_ITERATIVE_METHOD_HH
