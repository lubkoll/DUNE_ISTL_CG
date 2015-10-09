#ifndef DUNE_GENERIC_ITERATIVE_METHOD_HH
#define DUNE_GENERIC_ITERATIVE_METHOD_HH

#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include "voider.hh"
#include "util.hh"
#include "Mixins/maxSteps.hh"

namespace Dune
{
  /**
   * \cond
   */
  class InverseOperatorResult;

  inline std::ostream& operator<<(std::ostream& os, const InverseOperatorResult& res)
  {
    os << "InverseOperatorResult: converged         : " << std::boolalpha << res.converged << "\n";
    os << "                       iterations        : " << res.iterations << "\n";
    os << "                       reduction         : " << res.reduction << "\n";
    os << "                       convergence rate  : " << res.conv_rate << "\n";
    os << "                       elapsed time      : " << res.elapsed << " seconds\n";
    return os;
  }

  namespace GenericIterativeMethodDetail
  {
    template <class ToConnect, class Connector>
    using TryConnect = decltype(std::declval<Connector>().connect(std::declval<ToConnect>()));

    template <class Type>
    using TryMemFn_MinimalDecreaseAchieved = decltype( std::declval<Type>().minimalDecreaseAchieved() );

    template <class Type>
    using TryMemFn_Restart = decltype(std::declval<Type>().restart());

    template <class Type, class = void>
    struct HasMemFn_MinimalDecreaseAchieved : public std::false_type {};

    template <class Type>
    struct HasMemFn_MinimalDecreaseAchieved< Type , void_t< TryMemFn_MinimalDecreaseAchieved<Type> > > : public std::true_type {};

    template <class Type, bool hasMember = HasMemFn_MinimalDecreaseAchieved<Type>::value>
    struct BindMemFn_MinimalDecreaseAchieved
    {
      static std::function<bool()> apply(const Type& t)
      {
        return std::bind(&Type::minimalDecreaseAchieved,&t);
      }
    };

    template <class Type>
    struct BindMemFn_MinimalDecreaseAchieved<Type,false>
    {
      static std::function<bool()> apply(const Type&)
      {
        return std::function<bool()>{};
      }
    };

    template <class Type, class = void>
    struct OptionalRestart
    {
      static bool apply(const Type&)
      {
        return false;
      }
    };

    template <class Type>
    struct OptionalRestart< Type , void_t< TryMemFn_Restart<Type> > >
    {
      static bool apply(const Type& t)
      {
        return t.restart();
      }
    };

    template <class Type>
    bool optional_restart(const Type& t)
    {
      return OptionalRestart<Type>::apply(t);
    }


    template <class ToConnect, class Connector, class = void>
    struct ConnectIfPossible
    {
      static void apply(const ToConnect&, Connector&)
      {}
    };

    template <class ToConnect, class Connector>
    struct ConnectIfPossible< ToConnect , Connector , void_t<TryConnect<ToConnect,Connector> > >
    {
      static void apply(const ToConnect& toConnect, Connector& connector)
      {
        connector.connect(toConnect);
      }
    };


    template <class ToConnect, class Connector>
    void connect_if_possible(const ToConnect& toConnect, Connector& connector)
    {
      ConnectIfPossible<ToConnect,Connector>::apply(toConnect,connector);
    }


    template <class ToConnect, class Connector>
    void bind_connect_if_possible(const ToConnect& toConnect, Connector& connector)
    {
      auto relaxedTerminationCriterion = BindMemFn_MinimalDecreaseAchieved<ToConnect>::apply(toConnect);
      if( relaxedTerminationCriterion )
        connect_if_possible(relaxedTerminationCriterion,connector);
    }

    struct Empty{};

    template <class Step>
    using VerbosityIfNotDefined = access_t<std::conditional< std::is_base_of<Mixin::Verbosity,Step>::value , Empty , Mixin::Verbosity > >;
  }
  /**
   * \endcond
   */


  /*!
    @ingroup ISTL_Solvers
    @brief Generic wrapper for iterative methods.
   */
  template <class Step_,
            class TerminationCriterion_>
  class GenericIterativeMethod :
      public Step_ ,
      public InverseOperator<typename Step_::domain_type, typename Step_::range_type> ,
      public Mixin::MaxSteps ,
      public GenericIterativeMethodDetail::VerbosityIfNotDefined<Step_>
  {
  public:
    using Step = Step_;
    using TerminationCriterion = TerminationCriterion_;
    using domain_type = typename Step::domain_type;
    using range_type  = typename Step::range_type;
    using field_type  = field_t<domain_type>;
    using real_type   = real_t<Step>;

    /*!
      @param step object implementing one step of an iterative scheme
      @param terminate termination criterion
     */
    GenericIterativeMethod(Step step, TerminationCriterion terminate)
      : Step_(std::move(step)),
        terminate_(std::move(terminate))
    {
      terminate_.connect(*this);
      GenericIterativeMethodDetail::bind_connect_if_possible(terminate_,*this);
    }

    /*!
      @brief Apply loop solver to solve \f$Ax=b\f$.
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

        if( terminate_ ) break;
        if( GenericIterativeMethodDetail::optional_restart(*this) )
        {
          restoreInitialInput(x,b);
          Step::init(x,b);
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
      @brief Apply loop solver to solve \f$Ax=b\f$.
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
      @brief Apply loop solver to solve \f$Ax=b\f$.
      @param x initial iterate
      @param b initial right hand side
     */
    void apply(domain_type &x, range_type &b)
    {
      InverseOperatorResult res;
      apply(x,b,res);
    }

    /*!
      @brief Access termination criterion.
     */
    TerminationCriterion& terminationCriterion()
    {
      return terminate_;
    }

  private:
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