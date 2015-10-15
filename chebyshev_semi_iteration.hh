#ifndef DUNE_CHEBYSHEV_SEMI_ITERATION_HH
#define DUNE_CHEBYSHEV_SEMI_ITERATION_HH

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include "cg_solver.hh"
#include "residual_based_termination_criterion.hh"
#include "util.hh"

namespace Dune
{
  namespace ChebyshevSemiIterationSpec
  {
    //! Data object for the Chebyshev semi-iteration.
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
        if( !initialized_ )
          throw std::runtime_error("Uninitialized spectral bounds in chebyshev semi-iteration.");

        dx_ = std::make_unique<Domain>(x);
        *dx_ *= 0;
        b_ = std::make_unique<Range>(b);
        r_ = std::make_unique<Range>(b);
        A_->applyscaleadd(-1,x,*r_);
        Pr_ = std::make_unique<Domain>(x);
        P_->apply(*Pr_,*r_);
        sigma_ = sp_->dot(*r_,*Pr_);
        x_1_ = std::make_unique<Domain>(x); *x_1_ *= 0;
      }

      void reset(Domain& x, Range& b)
      {
        *dx_ *= 0;
        *Pr_ *= 0;
        step_ = 1;
        b_ = std::make_unique<Range>(b);
        r_ = std::make_unique<Range>(b);
        A_->applyscaleadd(-1,x,*r_);
        x_1_ = std::make_unique<Domain>(x); *x_1_ *= 0;
      }

      LinearOperator<Domain,Range>* A_;
      Preconditioner<Domain,Range>* P_;
      SeqScalarProduct<Domain> ssp_;
      ScalarProduct<Domain>* sp_;

      real_type spectralCenter_ = 0., spectralRadius_ = 0.;
      real_type alpha_ = 0., sigma_ = -1, beta_ = 0;
      std::unique_ptr<Domain> dx_ = nullptr, Pr_ = nullptr;
      std::unique_ptr<Range> b_ = nullptr, r_ = nullptr;
      std::unique_ptr<Domain> x_1_ = nullptr;
      unsigned step_ = 1;
      bool initialized_ = false;
    };


    //! Extends public interface of GenericStep for the Chebyshev semi-iteration.
    template <class Data>
    class Interface
    {
    public:
      template <class... Args>
      Interface(Args&&... args)
        : data_(std::forward<Args>(args)...)
      {}

      //! @brief Access norm of residual with respect to the norm induced by the preconditioner, i.e. \f$(r,Pr)\f$, where \f$r=b-Ax\f$.
      auto residualNorm() const
      {
        return data_.sigma_;
      }

      std::string name() const
      {
        return "Chebyshev Semi-Iteration";
      }

      template <class Real>
      void setSpectrum(Real center, Real radius)
      {
        data_.spectralCenter_ = center;
        data_.spectralRadius_ = radius;
        data_.initialized_ = true;
      }

      template <class Real>
      void setSpectralBounds(Real a, Real b)
      {
        setSpectrum( (a+b)/2 , (std::max(a,b) - std::min(a,b))/2 );
      }

      /*!
        \brief Sets spectral bounds for the case that \f$A\f$ is a mass matrix and a one-step Jacobi-preconditioner is used.

        In this case the spectrum of the preconditioned mass matrix is contained in \f$[0.5,2.5]\f$ (see @cite Wathen1987).

        @warning If you use a block-Jacobi instead of a Jacobi-preconditioner the spectral bounds are not correct any more. In this case increase the parameter halfSpectralDiameter.

        \param halfSpectralDiameter
       */
      void initializeForMassMatrix_TetrahedralQ1Elements(double halfSpectralDiameter = 1)
      {
        setSpectrum( 1 + 0.5*halfSpectralDiameter , halfSpectralDiameter );
      }

    protected:
      Data data_;
    };


    //! Compute parameters for the step computation in the Chebyshev semi-iteration.
    class ComputeStepParameters
    {
    public:
      template <class Data>
      void operator()(Data& data) const
      {
        data.sigma_ = data.sp_->dot(*data.Pr_,*data.r_);
        if( data.step_ == 1 )
        {
          data.beta_ = 0;
          data.alpha_ = -data.spectralCenter_;
          return;
        }

        if( data.step_ == 2 )
          data.beta_ = -0.5*data.spectralRadius_*data.spectralRadius_/data.spectralCenter_;
        else
          data.beta_ = 0.25*data.spectralRadius_*data.spectralRadius_/data.alpha_;

        data.alpha_ = -(data.spectralCenter_ + data.beta_);
      }
    };


    class UpdateIterate
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x) const
      {
        auto oldX = x;
        x *= data.spectralCenter_;
        x.axpy(1.,*data.Pr_);
        x.axpy(data.beta_,*data.x_1_);
        x *= -1/data.alpha_;
        *data.x_1_ = oldX;
      }
    };


    class ExplicitlyComputeResidual
    {
    public:
      template <class Data, class Domain>
      void operator()(Data& data, Domain& x) const
      {
        *data.r_ = *data.b_;
        data.A_->applyscaleadd(-1,x,*data.r_);
        data.P_->apply(*data.Pr_,*data.r_);
        data.sigma_ = data.sp_->dot(*data.r_,*data.Pr_);
      }
    };


    //! Step implementation for the chebshev semi-iteration.
    template <class Domain, class Range=Domain>
    using Step =
    GenericStep<Domain, Range,
      CGSpec::ApplyPreconditioner,
      GenericStepDetail::Ignore,
      ComputeStepParameters,
      GenericStepDetail::Ignore,
      UpdateIterate,
      ExplicitlyComputeResidual,
      Data<Domain,Range>,
      Interface
    >;
  }


  /**
   * @ingroup ISTL_Solvers
   * @brief Preconditioned chebyshev semi-iteration.
   *
   * Standard implementation based on a three-term recurrence and explicit computation of the residuals to avoid
   * accumulation of round of errors in the computation of the residuals. In contrast to Krylov methods this does not
   * slow down convergence (see @cite Gutknecht2002).
   *
   * @warning The Chebyshev semi-iteration requires bounds on the spectrum such as derived in @cite Wathen1987 for the case that a one-step Jacobi-preconditioner is used.
   */
  template <class Domain, class Range=Domain>
  using ChebyshevSemiIteration = GenericIterativeMethod< ChebyshevSemiIterationSpec::Step<Domain,Range> , KrylovTerminationCriterion::ResidualBased< real_t<Domain> > >;
}

#endif // DUNE_CHEBYSHEV_SEMI_ITERATION_HH
