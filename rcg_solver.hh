#ifndef DUNE_RCG_SOLVER_HH
#define DUNE_RCG_SOLVER_HH

#include <limits>
#include <utility>

#include "cg_solver.hh"
#include "generic_iterative_method.hh"
#include "generic_step.hh"
#include "operator_type.hh"
#include "relative_energy_termination_criterion.hh"
#include "tcg_solver.hh"
#include "util.hh"
#include "mixins/verbosity.hh"

namespace Dune
{
  namespace RCGSpec
  {
    //! Data object for the regularized conjugate gradient method.
    template <class Domain, class Range>
    struct Data : TCGSpec::Data<Domain,Range>
    {
      using real_type = real_t<Domain>;

      template <class... Args>
      Data(Args&&... args)
        : TCGSpec::Data<Domain,Range>(std::forward<Args>(args)...)
      {}

      void init(Domain& x, Range& b)
      {
        TCGSpec::Data<Domain,Range>::init(x,b);
        Pdx_ = std::make_unique<Range>(*this->r_);
      }

      void reset(Domain& x, Range& b)
      {
        TCGSpec::Data<Domain,Range>::reset(x,b);
        *Pdx_ = *this->r_;
        doRestart_ = false;
      }

      real_type theta_ = 0, dxPdx_ = 0;
      real_type minIncrease_ = 2, maxIncrease_ = 1000;
      std::unique_ptr<Range> Pdx_ = nullptr;
      bool doRestart_ = false;
    };


    //! Extends public interface of GenericStep for the regularized conjugate gradient method.
    template <class Data, class Name>
    class InterfaceImpl : public TCGSpec::InterfaceImpl<Data,Name>
    {
    public:
      template <class... Args>
      InterfaceImpl(Args&&... args)
        : TCGSpec::InterfaceImpl<Data,Name>(std::forward<Args>(args)...)
      {}

      //!* @brief Restart the regularized conjugate gradient method after regularization.
      bool restart() const
      {
        return data_.doRestart_;
      }

      /**
       * @brief Set minimal ratio for increasing the regularization parameter, i.e. \f$\frac{\theta_{new}}{\theta_{old}}>=minIncrease\f$.
       * @param minIncrease minimal ratio for increasing the regularization parameter
       */
      template <class Type>
      void setMinimalIncrease(Type minIncrease)
      {
        data_.minIncrease_ = minIncrease;
      }

      /**
       * @brief Set maximal ratio for increasing the regularization parameter, i.e. \f$\frac{\theta_{new}}{\theta_{old}}<=maxIncrease\f$.
       * @param maxIncrease maximal ratio for increasing the regularization parameter
       */
      template <class Type>
      void setMaximalIncrease(Type maxIncrease)
      {
        data_.maxIncrease_ = maxIncrease;
      }

    protected:
      using TCGSpec::InterfaceImpl<Data,Name>::data_;
    };

    /*! @cond */
    class Name
    {
    public:
      std::string name() const
      {
        return "Regularized Conjugate Gradients";
      }
    };
    /*! @endcond */

    //! Bind second template argument of RCG::InterfaceImpl to satisfy the interface of GenericStep.
    template <class Data>
    using Interface = InterfaceImpl<Data,Name>;


    //! Compute search direction.
    class SearchDirection : public CGSpec::SearchDirection
    {
    public:
      template <class Data>
      void operator()(Data& data) const
      {
        CGSpec::SearchDirection::operator()(data);

        // adjust energy norm of correction
        data.dxPdx_ = data.sp_->dot(*data.dx_,*data.Pdx_);
        data.dxAdx_ += data.theta_ * data.dxPdx_;
        // adjust preconditioned correction
        *data.Pdx_ *= data.beta_;
        *data.Pdx_ += *data.r_;
      }
    };


    //! Update residual.
    class UpdateResidual : public CGSpec::UpdateResidual
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x) const
      {
        CGSpec::UpdateResidual::operator()(data,x);
        data.r_->axpy(-data.alpha_*data.theta_,*data.Pdx_);
      }
    };


    //! Regularize if a direction of non-positive curvature is encountered.
    template <class real_type>
    class TreatNonconvexity :
        public Mixin::Eps<real_type>,
        public Mixin::Verbosity
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain&) const
      {
        if (data.dxAdx_ > 0 ) return;

        if( verbosityLevel() > 1 )
          std::cout << "    Regularizing at nonconvexity: " << data.dxAdx_ << std::endl;
        auto oldTheta = data.theta_ > 0 ? data.theta_ : this->eps();
        using std::abs;
        data.theta_ += (1-data.dxAdx_)/abs(data.dxPdx_);
        using std::min;
        using std::max;
        data.theta_ = min(max(data.minIncrease_*oldTheta,data.theta_),data.maxIncrease_*oldTheta);
        if( verbosityLevel() > 1 ) std::cout << "Updating regularization parameter from " << oldTheta << " to " << data.theta_ << std::endl;

        data.alpha_ = 0;
        data.operatorType_ = OperatorType::Indefinite;
        data.doRestart_ = true;
      }
    };


    //! Step implementation for the regularized conjugate gradient method.
    template <class Domain, class Range=Domain>
    using Step =
    GenericStep<Domain, Range,
      CGSpec::ApplyPreconditioner,
      SearchDirection,
      CGSpec::Scaling,
      TreatNonconvexity< real_t<Domain> >,
      CGSpec::UpdateIterate,
      UpdateResidual,
      Data<Domain,Range>,
      Interface
    >;
  }


  /*!
    @ingroup ISTL_Solvers
    @brief Regularized conjugate gradient method.

    Compute a descent direction for quadratic optimization problems of the
    form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a possibly indefinite linear operator.

    Regularizes if a conjugate search direction \f$\delta x\f$ of non-positive curvature is encountered, i.e. if \f$\delta xA\delta x\le 0\f$.
    In this case the operator \f$A\f$ is replaced by the operator \f$A+\theta P\f$, where \f$\theta\f$ is a monotone increasing regularization parameter.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::KrylovTerminationCriterion::ResidualBased or Dune::KrylovTerminationCriterion::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = KrylovTerminationCriterion::RelativeEnergyError>
  using RCGSolver = GenericIterativeMethod< RCGSpec::Step<Domain,Range> , TerminationCriterion< real_t<Domain> > >;
}

#endif // DUNE_RCG_SOLVER_HH
