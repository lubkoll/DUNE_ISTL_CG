#ifndef DUNE_TRCG_SOLVER_HH
#define DUNE_TRCG_SOLVER_HH

#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include <dune/common/typetraits.hh>

#include "cg_solver.hh"
#include "generic_iterative_method.hh"
#include "generic_step.hh"
#include "rcg_solver.hh"
#include "relative_energy_termination_criterion.hh"

namespace Dune
{
  namespace TRCGSpec
  {
    //! Data object for the truncated regularized conjugate gradient method.
    template <class Domain, class Range>
    struct Data : RCGSpec::Data<Domain,Range>
    {
      template <class... Args>
      Data(Args&&... args)
        : RCGSpec::Data<Domain,Range>(std::forward<Args>(args)...)
      {}

      std::function<bool()> minimalDecreaseAchieved_ = {};
    };

    /*! @cond */
    class Name
    {
    public:
      std::string name() const
      {
        return "Truncated Regularized Conjugate Gradients";
      }
    };
    /*! @endcond */


    //! Extends public interface of GenericStep for the truncated regularized conjugate gradient method.
    template <class Data>
    class Interface : public RCGSpec::InterfaceImpl<Data,Name>
    {
    public:
      template <class... Args>
      Interface(Args&&... args)
        : RCGSpec::InterfaceImpl<Data,Name>(std::forward<Args>(args)...)
      {}

      void connect(std::function<bool()> minimalDecreaseAchieved)
      {
        data_.minimalDecreaseAchieved_ = std::move(minimalDecreaseAchieved);
      }

      bool terminate() const
      {
        return data_.doTerminate_;
      }

    protected:
      using RCGSpec::InterfaceImpl<Data,Name>::data_;
    };


    //! Regularize or truncate at directions of non-positive curvature.
    template <class real_type>
    class TreatNonconvexity : public RCGSpec::TreatNonconvexity<real_type>
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x)
      {
        if( data.dxAdx_ > 0 ) return;

        assert(data.minimalDecreaseAchieved_);
        if( data.minimalDecreaseAchieved_() )
        {
          if( this->verbosityLevel() > 1 )
            std::cout << "    " << "Truncating at nonconvexity." << std::endl;
          data.alpha_ = 0;
          data.operatorType_ = OperatorType::Indefinite;
          data.doTerminate_ = true;
          return;
        }

        RCGSpec::TreatNonconvexity<real_type>::operator()(data,x);
      }
    };

    //! Step implementation for the truncated regularized conjugate gradient method.
    template <class Domain, class Range=Domain>
    using Step =
    GenericStep<Domain, Range,
      CGSpec::ApplyPreconditioner,
      RCGSpec::SearchDirection,
      CGSpec::Scaling,
      TreatNonconvexity< real_t<Domain> >,
      CGSpec::UpdateIterate,
      RCGSpec::UpdateResidual,
      Data<Domain,Range>,
      Interface
    >;
  }


  /*!
    @ingroup ISTL_Solvers
    @brief Truncated regularized conjugate gradient method (see @cite Lubkoll2015a).

    Compute a descent direction for quadratic optimization problems of the
    form \f$ q(x)=\frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a possibly indefinite linear operator.

    Combines the regularization strategy of RCG with TCG if a conjugate search direction \f$\delta x\f$ of non-positive curvature is encountered, i.e. if \f$\delta xA\delta x\le 0\f$.

    Suppose that \f$q\f$ is a local model of a nonconvex optimization problem. Far from the solution seach directions do not need to be computed overly accurate, say up to some relative
    accuracy \f$\delta_{min}\f$. Close to the solution we need to compute more accurate corrections, say up to some relative accuracy \f$\delta\f$. In this setting the TRCG method
    treats directions of non-positive curvature as follows:
     - If \f$\delta < \delta_{min}\f$ and the current iterate is acceptable with respect to \f$\delta_{min}\f$ then we conclude that we are still far from the solution, (where \f$A\f$ is positive definite)
       and accept the computed iterate.
     - Else the regularization strategy of RCG is applied.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::KrylovTerminationCriterion::RelativeEnergyError, must provide a member function bool minimalDecreaseAchieved())
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = KrylovTerminationCriterion::RelativeEnergyError>
  using TRCGSolver = GenericIterativeMethod< TRCGSpec::Step<Domain,Range> , TerminationCriterion< real_t<Domain> > >;
}

#endif // DUNE_TRCG_SOLVER_HH
