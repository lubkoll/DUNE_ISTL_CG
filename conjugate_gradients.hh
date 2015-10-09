#ifndef DUNE_CONJUGATE_GRADIENTS_HH
#define DUNE_CONJUGATE_GRADIENTS_HH

#include <dune/common/static_assert.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include "conjugate_gradient_step.hh"
#include "generic_iterative_method.hh"
#include "relative_energy_termination_criterion.hh"
#include "util.hh"

namespace Dune
{
  /*!
    @ingroup ISTL_Solvers
    @brief Conjugate gradient method (see @cite Hestenes1952).

    Solves quadratic optimization problems of the form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a positive definite linear operator.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::Termination::ResidualBased or Dune::Termination::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = Termination::RelativeEnergyError >
  using CG   = GenericIterativeMethod< CGStepImpl<Domain,Range,CGDetail::CGBase>    , TerminationCriterion< real_t<Domain> > >;

  /*!
    @ingroup ISTL_Solvers
    @brief Regularized conjugate gradient method.

    Compute a descent direction for quadratic optimization problems of the
    form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a possibly indefinite linear operator.

    Regularizes if a conjugate search direction \f$\delta x\f$ of non-positive curvature is encountered, i.e. if \f$\delta xA\delta x\le 0\f$.
    In this case the operator \f$A\f$ is replaced by the operator \f$A+\theta P\f$, where \f$\theta\f$ is a monotone increasing regularization parameter.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::Termination::ResidualBased or Dune::Termination::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = Termination::RelativeEnergyError >
  using RCG  = GenericIterativeMethod< CGStepImpl<Domain,Range,CGDetail::RCGBase> , TerminationCriterion< real_t<Domain> > >;

  /*!
    @ingroup ISTL_Solvers
    @brief Truncated conjugate gradient method.

    Compute a descent direction for quadratic optimization problems of the
    form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a possibly indefinite linear operator.

    Terminates if a conjugate search direction \f$\delta x\f$ of non-positive curvature is encountered, i.e. if \f$\delta xA\delta x\le 0\f$.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::Termination::ResidualBased or Dune::Termination::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = Termination::RelativeEnergyError >
  using TCG  = GenericIterativeMethod< CGStepImpl<Domain,Range,CGDetail::TCGBase>  , TerminationCriterion< real_t<Domain> > >;

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
    @tparam TerminationCriterion termination criterion (such as Dune::Termination::ResidualBased or Dune::Termination::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = Termination::RelativeEnergyError >
  using TRCG = GenericIterativeMethod< CGStepImpl<Domain,Range,CGDetail::TRCGBase> , TerminationCriterion< real_t<Domain> > >;


  /*!
    @ingroup ISTL_Solvers
    @brief Generate conjugate gradient method.

    Solves equations of the form \f$PAx=Pb\f$, where \f$A:\ X\mapsto Y\f$ is a linear operator and
    \f$ P:\ Y\mapsto X\f$ a preconditioner.

    Usage:
    @code{.cpp}
    auto cg   = make_cg<Dune::CG  ,Dune::Termination::ResidualBased>(A,P,sp,...);
    auto rcg  = make_cg<Dune::RCG ,Dune::Termination::ResidualBased>(A,P,sp,...);
    auto tcg  = make_cg<Dune::TCG ,Dune::Termination::ResidualBased>(A,P,sp,...);
    auto trcg = make_cg<Dune::TRCG,Dune::Termination::ResidualBased>(A,P,sp,...);
    @endcond

    @param A linear operator
    @param P preconditioner
    @param sp scalar product
    @param accuracy relative accuracy
    @param nSteps maximal number of steps
    @param verbosityLevel =1: print final statistics, =2: print information in each iteration
    @param eps maximal attainable accuracy
    @tparam CGType conjugate gradient variant (=CG,RCG,TCG or TRCG)
    @tparam TerminationCriterion termination criterion (such as Dune::Termination::ResidualBased or Dune::Termination::RelativeEnergyError)
    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
   */
  template <template <class,class,template <class> class> class CGType,
            template <class> class TerminationCriterion,
            class Domain, class Range, class real_type = real_t<Domain> >
  auto make_cg(LinearOperator<Domain,Range>& A,
               Preconditioner<Domain,Range>& P,
               ScalarProduct<Domain>& sp,
               real_type accuracy = 1e-15, unsigned nSteps = 1000,
               unsigned verbosityLevel = 0, real_type eps = 1e-15)
  {
    using MyCG = CGType<Domain,Range,TerminationCriterion>;
    TerminationCriterion< real_t<Domain> > terminationCriterion;
    terminationCriterion.setRelativeAccuracy(accuracy);
    terminationCriterion.setEps(eps);
    auto cg = MyCG{ typename MyCG::Step{ A,P,sp } , std::move(terminationCriterion) };
    cg.setMaxSteps(nSteps);
    cg.setVerbosityLevel(verbosityLevel);
    return cg;
  }
}

#endif // DUNE_CONJUGATE_GRADIENTS_HH
