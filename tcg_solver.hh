#ifndef DUNE_TCG_SOLVER_HH
#define DUNE_TCG_SOLVER_HH

#include <iostream>
#include <utility>

#include "cg_solver.hh"
#include "generic_iterative_method.hh"
#include "generic_step.hh"
#include "operator_type.hh"
#include "relative_energy_termination_criterion.hh"
#include "util.hh"
#include "mixins/verbosity.hh"

namespace Dune
{
  namespace TCGSpec
  {
    //! Data object for the truncated conjugate gradient method.
    template <class Domain, class Range>
    struct Data : CGSpec::Data<Domain,Range>
    {
      template <class... Args>
      Data(Args&&... args)
        : CGSpec::Data<Domain,Range>(std::forward<Args>(args)...)
      {}

      void reset(Domain& x, Domain& b)
      {
        CGSpec::Data<Domain,Range>::reset(x,b);
        doTerminate_ = false;
      }

      OperatorType operatorType_ = OperatorType::PositiveDefinite;
      bool doTerminate_ = false;
      bool performBlindUpdate_ = true;
    };

    /*! @cond */
    class Name
    {
    public:
      std::string name() const
      {
        return "Truncated Conjugate Gradients";
      }
    };
    /*! @endcond */

    //! Extends public interface of GenericStep for the truncated conjugate gradient method.
    template <class Data, class Name>
    class InterfaceImpl : public CGSpec::InterfaceImpl<Data,Name>
    {
    public:
      template <class... Args>
      InterfaceImpl(Args&&... args)
        : CGSpec::InterfaceImpl<Data,Name>(std::forward<Args>(args)...)
      {}

      bool terminate() const
      {
        return data_.doTerminate_;
      }

      bool operatorIsPositiveDefinite() const
      {
        return data_.operatorType_ == OperatorType::PositiveDefinite;
      }

      void setPerformBlindUpdate(bool blindUpdate = true)
      {
        data_.performBlindUpdate_ = blindUpdate;
      }

    protected:
      using CGSpec::InterfaceImpl<Data,Name>::data_;
    };


    //! Bind second template argument of TCG::InterfaceImpl to satisfy the interface of GenericStep.
    template <class Data>
    using Interface = InterfaceImpl<Data,Name>;


    //! Truncate at directions of non-positive curvature.
    class TreatNonconvexity : public Mixin::Verbosity
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x) const
      {
        if( data.dxAdx_ > 0 ) return;

        if( verbosityLevel() > 1 )
          std::cout << "    " << "Truncating at nonconvexity" << std::endl;
        // At least do something to retain a little chance to get out of the nonconvexity. If a nonconvexity is encountered in the first step something probably went wrong
        // elsewhere. Chances that a way out of the nonconvexity can be found are small in this case.
        if( data.performBlindUpdate_ )
          x += *data.dx_;
        data.doTerminate_ = true;
        data.operatorType_ = OperatorType::Indefinite;
      }
    };


    //! Step implementation for the truncated conjugate gradient method.
    template <class Domain, class Range=Domain>
    using Step =
    GenericStep<Domain, Range,
      CGSpec::ApplyPreconditioner,
      CGSpec::SearchDirection,
      CGSpec::Scaling,
      TreatNonconvexity,
      CGSpec::UpdateIterate,
      CGSpec::UpdateResidual,
      Data<Domain,Range>,
      Interface
    >;
  }


  /*!
    @ingroup ISTL_Solvers
    @brief Truncated conjugate gradient method.

    Compute a descent direction for quadratic optimization problems of the
    form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a possibly indefinite linear operator.

    Terminates if a conjugate search direction \f$\delta x\f$ of non-positive curvature is encountered, i.e. if \f$\delta xA\delta x\le 0\f$.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::KrylovTerminationCriterion::ResidualBased or Dune::KrylovTerminationCriterion::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = KrylovTerminationCriterion::RelativeEnergyError>
  using TCGSolver = GenericIterativeMethod< TCGSpec::Step<Domain,Range> , TerminationCriterion< real_t<Domain> > >;
}

#endif // DUNE_TCG_SOLVER_HH
