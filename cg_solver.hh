#ifndef DUNE_CG_HH
#define DUNE_CG_HH

#include <memory>
#include <utility>

#include "generic_iterative_method.hh"
#include "generic_step.hh"
#include "relative_energy_termination_criterion.hh"
#include "mixins/iterativeRefinements.hh"

namespace Dune
{
  namespace CGSpec
  {
    //! Data object for the conjugate gradient method.
    template <class Domain, class Range>
    struct Data
    {
      using real_type = real_t<Domain>;

      template <class LinOp, class Prec, class SP>
      Data(LinOp& A, Prec& P, SP& sp)
        : A_(&A), P_(&P), ssp_(), sp_(&sp)
      {
        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
      }

      template <class LinOp, class Prec>
      Data(LinOp& A,  Prec& P)
        : A_(&A), P_(&P), ssp_(), sp_(&ssp_)
      {
        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
      }

      Data(Data&&) = default;
      Data& operator=(Data&&) = default;

      Data(const Data& data)
        : A_(data.A_), P_(data.P_), ssp_(), sp_(data.sp_)
      {}

      Data& operator=(const Data& data)
      {
        A_ = data.A_;
        P_ = data.P_;
        sp_ = data.sp_;
      }

      void init(Domain& x, Range& b)
      {
        dx_ = nullptr;
        Adx_ = std::make_unique<Range>(b); *Adx_ *= 0.,
        Pr_ = std::make_unique<Domain>(x); *Pr_ *= 0.;
        r_ = &b;//std::make_unique<Range>(b);
        A_->applyscaleadd(-1.,x,*r_);
      }

      void reset(Domain& x, Range& b)
      {
        dx_ = nullptr;
        *Adx_ *= 0.;
        *Pr_ *= 0.;
        r_ = &b;//std::make_unique<Range>(b);
        A_->applyscaleadd(-1.,x,*r_);
      }

      LinearOperator<Domain,Range>* A_;
      Preconditioner<Domain,Range>* P_;
      SeqScalarProduct<Domain> ssp_;
      ScalarProduct<Domain>* sp_;

      std::unique_ptr<Range> Adx_ = nullptr;
      Range* r_ = nullptr;
      std::unique_ptr<Domain> Pr_ = nullptr, dx_ = nullptr;
      real_type alpha_ = -1, beta_ = -1, sigma_ = -1, dxAdx_ = -1;
    };

    //! @cond
    class Name
    {
    public:
      std::string name() const
      {
        return "Conjugate Gradients";
      }
    };
    //! @endcond


    //! Extends public interface of GenericStep for the conjugate gradient method.
    template <class Data, class NameOfAlgorithm = Name>
    class InterfaceImpl : public NameOfAlgorithm
    {
    public:
      template <class... Args>
      InterfaceImpl(Args&&... args)
        : data_(std::forward<Args>(args)...)
      {}

      //! @brief Access scaling for the conjugate search direction, i.e. \f$\frac{(r,Pr)}{(\delta x,A\delta x)}\f$
      auto alpha() const
      {
        return data_.alpha_;
      }

      //! @brief Access length of conjugate search direction with respect to the energy norm, i.e. \f$(\delta x,A\delta x)\f$.
      auto length() const
      {
        return data_.dxAdx_;
      }

      //! @brief Access norm of residual with respect to the norm induced by the preconditioner, i.e. \f$(r,Pr)\f$, where \f$r=b-Ax\f$.
      auto preconditionedResidualNorm() const
      {
        return data_.sigma_;
      }

      //! @brief Access norm of residual with respect to the employed scalar product (in general the l2-norm), i.e. \f$\|r\|_2\f$, where \f$r=b-Ax\f$.
      auto residualNorm()
      {
        return data_.sp_->norm(*data_.r_);
      }

    protected:
      Data data_;
    };

    //! Bind second template argument of CG::InterfaceImpl to satisfy the interface of GenericStep.
    template <class Data>
    using Interface = InterfaceImpl<Data,Name>;


    //! Apply preconditioner, possibly with iterative refinements.
    class ApplyPreconditioner
        : public Mixin::IterativeRefinements
    {
    public:
      template <class Data>
      void operator()(Data& data) const
      {
        assert(data.r_ != nullptr);
        assert(data.Pr_ != nullptr);

        data.P_->apply(*data.Pr_,*data.r_);

        if( iterativeRefinements() > 0 )
        {
          auto r2 = *data.r_;
          auto dQr = *data.Pr_;
          for(auto i=0u; i<iterativeRefinements(); ++i)
          {
            data.A_->applyscaleadd(-1.,*data.Pr_,r2);
            data.P_->apply(dQr,r2);
            *data.Pr_ += dQr;
          }
        }

        using std::abs;
        if( data.sigma_ < 0 )
          data.sigma_ = abs( data.sp_->dot(*data.r_,*data.Pr_) );
      }

      template <class Data, class Domain, class Range>
      void pre(Data& data, Domain& x, Range& b) const
      {
        data.P_->pre(x,b);
      }

      template <class Data, class Domain>
      void post(Data& data, Domain& x) const
      {
        data.P_->post(x);
      }
    };


    //! Compute search direction for the conjugate gradient method.
    class SearchDirection
    {
    public:
      template <class Data>
      void operator()(Data& data) const
      {
        if( data.dx_ == nullptr)
        {
          data.dx_ = std::make_unique< std::decay_t<decltype(*data.Pr_)> >(*data.Pr_);
          computeInducedStepLength(data);
          return;
        }

        using std::abs;
        auto newSigma = abs( data.sp_->dot(*data.r_,*data.Pr_) );
        data.beta_ = newSigma/data.sigma_;
        *data.dx_ *= data.beta_; *data.dx_ += *data.Pr_;
        data.sigma_ = newSigma;

        computeInducedStepLength(data);
      }

    private:
      template <class Data>
      void computeInducedStepLength(Data& data) const
      {
        data.A_->apply(*data.dx_,*data.Adx_);
        data.dxAdx_ = data.sp_->dot(*data.dx_,*data.Adx_);
      }
    };


    //! Compute scaling of the search direction for the conjugate gradient method.
    class Scaling
    {
    public:
      template <class Data>
      void operator()(Data& data) const
      {
        data.alpha_ = data.sigma_/data.dxAdx_;
      }
    };


    class UpdateIterate
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x) const
      {
        x.axpy(data.alpha_,*data.dx_);
      }
    };


    //! Adjust residual.
    class UpdateResidual
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain&) const
      {
        data.r_->axpy(-data.alpha_,*data.Adx_);
      }
    };


    //! Step implementation for the conjugate gradient method.
    template <class Domain, class Range=Domain>
    using Step =
    GenericStep<Domain, Range,
      ApplyPreconditioner,
      SearchDirection,
      Scaling,
      GenericStepDetail::Ignore,
      UpdateIterate,
      UpdateResidual,
      Data<Domain,Range>,
      Interface
    >;
  }


  /*!
    @ingroup ISTL_Solvers
    @brief Conjugate gradient method (see @cite Hestenes1952).

    Solves quadratic optimization problems of the form \f$ \frac{1}{2}x^T Ax - b^T x \f$, where \f$A:\ X\mapsto Y\f$ is a positive definite linear operator.

    @tparam Domain domain space \f$X\f$
    @tparam Range range space \f$Y\f$
    @tparam TerminationCriterion termination criterion (such as Dune::KrylovTerminationCriterion::ResidualBased or Dune::KrylovTerminationCriterion::RelativeEnergyError (default))
   */
  template <class Domain, class Range,
            template <class> class TerminationCriterion = KrylovTerminationCriterion::RelativeEnergyError>
  using MyCGSolver = GenericIterativeMethod< CGSpec::Step<Domain,Range> , TerminationCriterion< real_t<Domain> > >;


  /*!
    @ingroup ISTL_Solvers
    @brief Generate conjugate gradient method.

    Solves equations of the form \f$PAx=Pb\f$, where \f$A:\ X\mapsto Y\f$ is a linear operator and
    \f$ P:\ Y\mapsto X\f$ a preconditioner.

    Usage:
    @code{.cpp}
    auto cg   = make_cg<Dune::CG  ,Dune::KrylovTerminationCriterion::ResidualBased>(A,P,sp,...);
    auto rcg  = make_cg<Dune::RCG ,Dune::KrylovTerminationCriterion::ResidualBased>(A,P,sp,...);
    auto tcg  = make_cg<Dune::TCG ,Dune::KrylovTerminationCriterion::ResidualBased>(A,P,sp,...);
    auto trcg = make_cg<Dune::TRCG,Dune::KrylovTerminationCriterion::ResidualBased>(A,P,sp,...);
    @endcond

    @param A linear operator
    @param P preconditioner
    @param sp scalar product
    @param accuracy relative accuracy
    @param nSteps maximal number of steps
    @param verbosityLevel =1: print final statistics, =2: print information in each iteration
    @param eps maximal attainable accuracy
    @tparam CGType conjugate gradient variant (=CG,RCG,TCG or TRCG)
    @tparam TerminationCriterion termination criterion (such as Dune::KrylovTerminationCriterion::ResidualBased or Dune::KrylovTerminationCriterion::RelativeEnergyError)
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

#endif // DUNE_CG_HH
