#ifndef DUNE_CHEBYSHEV_SEMI_ITERATION_HH
#define DUNE_CHEBYSHEV_SEMI_ITERATION_HH

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include <dune/common/typetraits.hh>
#include "cg_solver.hh"
#include "residual_based_termination_criterion.hh"

namespace Dune
{
  /**
   * @brief One step of the Chebyshev semi-iteration.
   *
   * @note Requires spectral bounds to be provided by one of the methods setSpectum(), setSpectralBounds() or initializeForMassMatrix_TetrahedralQ1Elements().
   */
  template < class Domain, class Range = Domain >
  class ChebyshevSemiIterationStep
  {
  public:
    //! type of the domain space
    using domain_type = Domain;
    //! type of the range space
    using range_type = Range;
    //! underlying field type
    using field_type = field_t<Domain>;
    //! corresponding real type (same as real type for real spaces, differs for complex spaces)
    using real_type = real_t<Domain>;

    /**
     * @brief Constructor.
     *
     * @param A linear operator
     * @param P preconditioner (often a Jacobi preconditioner)
     * @param sp scalar product
     */
    template <class LinOp, class Prec, class SP>
    ChebyshevSemiIterationStep(LinOp& A, Prec& P, SP& sp)
      : A_(A), P_(P), ssp_(), sp_(sp)
    {
//        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
//        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
    }

    /**
     * @brief Constructor.
     *
     * @param A linear operator
     * @param P preconditioner (often a Jacobi preconditioner)
     */
    template <class LinOp, class Prec>
    ChebyshevSemiIterationStep(LinOp& A,  Prec& P)
      : A_(A), P_(P), ssp_(), sp_(ssp_)
    {
//        static_assert( LinOp::category == Prec::category , "Linear operator and preconditioner are required to belong to the same category!" );
//        static_assert( LinOp::category == SolverCategory::sequential , "Linear operator must be sequential!" );
    }

    //! Initialization phase, initialize storage and apply preprocessing phase of the preconditioner.
    void init( Domain& x, Range& b )
    {
      P_.pre( x, b );
      reset( x, b );
    }

    //! Postprocessing phase, apply postprocessing of the preconditioner.
    void postProcess( Domain& x )
    {
      P_.post( x );
    }

    //! Reset internal storage.
    void reset( Domain& x, Range& b )
    {
      if( !initialized_ )
        throw std::runtime_error("Uninitialized spectral bounds in chebyshev semi-iteration.");

      if( dx_ == nullptr)
      {
        dx_ = std::unique_ptr<Domain>( new Domain(x) );
        Pr_ = std::unique_ptr<Domain>( new Domain(x) );
        x1_ = std::unique_ptr<Domain>( new Domain(x) );
        r_ = std::unique_ptr<Range>( new Range(b) );
      }

      A_.applyscaleadd(-1,x,*r_);
      P_.apply(*Pr_,*r_);

      sigma_ = sp_.dot(*r_,*Pr_);
      *x1_ *= 0;
      *dx_ *= 0;
      step_ = 1;
    }

    //! Perform one step of the Chebyshev semi-iteration.
    void compute( Domain& x, Range& b )
    {
      // compute parameters
      if( step_ == 1 )
      {
        beta_ = 0;
        alpha_ = -spectralCenter_;
      }
      else
      {
        if( step_ == 2 )
          beta_ = -0.5*spectralRadius_*spectralRadius_/spectralCenter_;
        else
          beta_ = 0.25*spectralRadius_*spectralRadius_/alpha_;

        alpha_ = -( spectralCenter_ + beta_ );
      }

      // update iterate
      auto oldX = x;
      x *= spectralCenter_;
      x.axpy( 1, *Pr_ );
      x.axpy( beta_, *x1_ );
      x *= -1/alpha_;
      *x1_ = oldX;

      // compute residual
      *r_ = b;
      A_.applyscaleadd( -1, x, *r_ );

      // apply preconditioner
      P_.apply( *Pr_, *r_ );
      sigma_ = sp_.dot( *r_, *Pr_ );
    }

    std::string name() const
    {
      return "Chebyshev Semi-Iteration";
    }

    /**
     * @brief Provide information about the spectrum.
     * @param center center of the spectrum
     * @param radius radius of the spectrum
     */
    void setSpectrum(real_t<Domain> center, real_t<Domain> radius)
    {
      spectralCenter_ = center;
      spectralRadius_ = radius;
      initialized_ = true;
    }

    /**
     * @brief Provide information about the spectrum.
     * @param a lower, resp. upper bound of the spectrum
     * @param b upper, resp. lower bound of the spectrum
     */
    void setSpectralBounds(real_t<Domain> a, real_t<Domain> b)
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
      setSpectrum( 0.5 + halfSpectralDiameter , halfSpectralDiameter );
    }

    //! @brief Access norm of residual with respect to the norm induced by the preconditioner, i.e. \f$(r,Pr)\f$, where \f$r=b-Ax\f$.
    double residualNorm() const
    {
      return sqrt(sigma_);
    }

  private:
    real_t<Domain> spectralCenter_ = 0., spectralRadius_ = 0.;
    real_t<Domain> alpha_ = 0., sigma_ = -1, beta_ = 0;
    std::unique_ptr<Domain> dx_ = nullptr, Pr_ = nullptr, x1_ = nullptr;
    std::unique_ptr<Range> r_ = nullptr;
    unsigned step_ = 1;
    bool initialized_ = false;

    LinearOperator<Domain,Range>& A_;
    Preconditioner<Domain,Range>& P_;
    SeqScalarProduct<Domain> ssp_;
    ScalarProduct<Domain>& sp_;
  };


  /**
   * @ingroup ISTL_Solvers
   * @brief Preconditioned chebyshev semi-iteration.
   *
   * Standard implementation based on a three-term recurrence and explicit computation of the residuals to avoid
   * accumulation of round of errors.
   *
   * If spectral bounds are available, then the Chebyshev semi-iteration, with a fixed size of steps, provides a linear preconditioner (see @cite Gutknecht2002).
   *
   * @warning The Chebyshev semi-iteration requires bounds on the spectrum such as derived in @cite Wathen1987 for the case that a one-step Jacobi-preconditioner is used.
   */
  template <class Domain, class Range=Domain>
  using ChebyshevSemiIteration = GenericIterativeMethod< ChebyshevSemiIterationStep<Domain,Range>, KrylovTerminationCriterion::ResidualBased< real_t<Domain> > >;
}

#endif // DUNE_CHEBYSHEV_SEMI_ITERATION_HH
