#ifndef DUNE_CONJUGATE_GRADIENT_STEP_HH
#define DUNE_CONJUGATE_GRADIENT_STEP_HH

#include <functional>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>

#include "Mixins/eps.hh"
#include "Mixins/iterativeRefinements.hh"
#include "Mixins/verbosity.hh"
#include "util.hh"

namespace Dune
{
  /*! @cond */
  namespace CGDetail
  {
    enum class OperatorType { PositiveDefinite, Indefinite };


    template <class domain_type,class range_type>
    class CGBase :
        public Mixin::Verbosity ,
        public Mixin::Eps< real_t<range_type> >
    {
      using real_type = real_t<range_type>;
    public:
      CGBase() = default;

      CGBase(ScalarProduct<domain_type>&)
      {}

      bool isPositiveDefinite() const
      {
        return operatorType_ == OperatorType::PositiveDefinite;
      }

    protected:
      void init(const range_type&)
      {}

      void reset()
      {}

      void regularize(real_type&, const domain_type&)
      {}

      void adjustRegularizedResidual(real_type , range_type&) const
      {}

      /*!
        @brief Handle search directions of non-positive curvature, i.e. directions \f$q\f$ with \f$qAq <= 0\f$.
        @throws if a direction of non-positive curvature is encountered in the standard conjugate gradient method.
       */
      void treatNonconvexity(real_type qAq, domain_type&, const domain_type&) const
      {
        if( verbose() )
        {
          std::cout << "    " << "CG: Direction of non-positive curvature encountered in standard CG Implementation!" << std::endl;
          std::cout << "    " << "CG: Either something is wrong with your operator or you should use TCG, RCG or HCG. Terminating CG!" << std::endl;
        }

        throw std::runtime_error("Non-positive curvature encountered in conjugate gradient method.");
      }

      void adjustPreconditionedSearchDirection(real_type, const range_type&)
      {}

      std::string name() const
      {
        return "Conjugate Gradients";
      }

      void setOperatorType(OperatorType operatorType)
      {
        operatorType_ = operatorType;
      }

      OperatorType operatorType_ = OperatorType::PositiveDefinite;
    };

    template <class domain_type,class range_type>
    class TCGBase : public CGBase<domain_type,range_type>
    {
      using Base = CGBase<domain_type,range_type>;
      using real_type = real_t<range_type>;
    public:
      TCGBase() = default;
      TCGBase(ScalarProduct<domain_type>&)
      {}

    protected:
      using Base::operatorType_;
      void init(const range_type&)
      {}

      void reset()
      {
        operatorType_ = OperatorType::PositiveDefinite;
        firstIteration_ = true;
      }

      void regularize(real_type&,const domain_type&)
      {
        firstIteration_ = false;
      }

      /*!
        @brief Handle search directions of negative curvature, i.e. directions \f$q\f$ with \f$qAq < 0\f$.
       */
      void treatNonconvexity(real_type, domain_type& x, const domain_type& q)
      {
        // At least do something to retain a little chance to get out of the nonconvexity. If a nonconvexity is encountered in the first step something probably went wrong
        // elsewhere. Chances that a way out of the nonconvexity can be found are small in this case.
        if( firstIteration_ ) x += q;
        this->setOperatorType(OperatorType::Indefinite);
      }

      std::string name() const
      {
        return "Truncated Conjugate Gradients";
      }

      bool firstIteration_ = true;
    };


    template <class domain_type,class range_type>
    class RCGBase :
        public CGBase<domain_type,range_type>
    {
      using Base = CGBase<domain_type,range_type>;
      using real_type = real_t<range_type>;
    public:
      RCGBase()
        : ssp_() , sp_(ssp_)
      {}

      RCGBase(ScalarProduct<domain_type>& sp)
        : sp_(sp)
      {}


      /**
       * @brief Restart the regularized conjugate gradient method after regularization.
       */
      bool restart() const
      {
        return operatorType_ == OperatorType::Indefinite;
      }

      /**
       * @brief Set minimal ratio for increasing the regularization parameter, i.e. \f$\frac{\theta_{new}}{\theta_{old}}>=minIncrease\f$.
       * @param minIncrease minimal ratio for increasing the regularization parameter
       */
      void setMinimalIncrease(real_type minIncrease)
      {
        minIncrease_ = minIncrease;
      }

      /**
       * @brief Set maximal ratio for increasing the regularization parameter, i.e. \f$\frac{\theta_{new}}{\theta_{old}}<=maxIncrease\f$.
       * @param maxIncrease maximal ratio for increasing the regularization parameter
       */
      void setMaximalIncrease(real_type maxIncrease)
      {
        maxIncrease_ = maxIncrease;
      }

    protected:
      using Base::operatorType_;

      void init(const range_type& residual)
      {
        Pdx_ = std::make_unique<range_type>(residual);
      }

      /**
       * @brief Reset regularization parameter, i.e. set \f$\theta=0\f$.
       */
      void reset()
      {
        theta_ = 0;
        this->setOperatorType(OperatorType::PositiveDefinite);
      }

      /**
       * @brief Replace \f$qAq\f$ with \f$qAq+\theta*qPq\f$, where \f$\theta\f$ is the regularization parameter.
       * @param qAq length of search direction with respect to A
       * @param q conjugate search direction
       * @param sp scalar product
       */
      void regularize(real_type& qAq, const domain_type& q)
      {
        dxPdx_ = sp_.dot(q,*Pdx_);
        qAq += theta_ * dxPdx_;
      }

      /**
       * @brief Replace \f$r\f$ with the regularized residual \f$r-\alpha*\theta*Pq\f$, where \f$\theta\f$ is the regularization parameter.
       * @param alpha step length parameter from the conjugate gradient method
       * @param Pq preconditioned search direction
       * @param r residual
       */
      void adjustRegularizedResidual(real_type alpha, range_type& r) const
      {
        r.axpy(-alpha*theta_,*Pdx_);
      }

      void treatNonconvexity(real_type qAq, domain_type&, const domain_type&)
      {
        updateRegularizationParameter(qAq);
        this->setOperatorType(OperatorType::Indefinite);
      }

      void adjustPreconditionedSearchDirection(real_type beta, const range_type& residual)
      {
        *Pdx_ *= beta; *Pdx_ += residual;
      }

      std::string name() const
      {
        return "Regularized Conjugate Gradients";
      }

    private:
      void updateRegularizationParameter(real_type qAq)
      {
        real_type oldTheta = theta_ > 0 ? theta_ : real_type(this->eps());
        using std::abs;
        theta_ += (1-qAq)/abs(dxPdx_);
        if( this->verbosityLevel() > 1 ) std::cout << "Computed regularization parameter: " << theta_ << std::endl;
        using std::min;
        using std::max;
        theta_ = min(max(minIncrease_*oldTheta,theta_),maxIncrease_*oldTheta);
        if( this->verbosityLevel() > 1 ) std::cout << "Updating regularization parameter from " << oldTheta << " to " << theta_ << std::endl;
      }

      SeqScalarProduct<domain_type> ssp_;
      ScalarProduct<domain_type>& sp_;
      real_type theta_ = 0, dxPdx_ = 0;
      real_type minIncrease_ = 2, maxIncrease_ = 1000;
      std::unique_ptr<range_type> Pdx_ = nullptr;
    };

    template <class domain_type,class range_type>
    class TRCGBase : public RCGBase<domain_type,range_type>
    {
      using RCGBase<domain_type,range_type>::operatorType_;
      using real_type = real_t<range_type>;
    public:
      TRCGBase() = default;

      TRCGBase(ScalarProduct<domain_type>& sp)
        : RCGBase<domain_type,range_type>(sp)
      {}

      void connect(std::function<bool()> minimalDecreaseAchieved)
      {
        minimalDecreaseAchieved_ = minimalDecreaseAchieved;
      }

    protected:
      void treatNonconvexity(real_type qAq, domain_type& x, const domain_type& q)
      {
        if( minimalDecreaseAchieved_() )
        {
          // At least do something to retain a little chance to get out of the nonconvexity. If a nonconvexity is encountered in the first step something probably went wrong
          // elsewhere. Chances that a way out of the nonconvexity can be found are small in this case.
          this->setOperatorType(OperatorType::Indefinite);
          return;
        }

        RCGBase<domain_type,range_type>::treatNonconvexity(qAq,x,q);
      }

      std::string name() const
      {
        return "Truncated Regularized Conjugate Gradients";
      }

    private:
      real_type theta_ = 0.;
      std::function<bool()> minimalDecreaseAchieved_;
    };
  }
  /*! @endcond */


  /*!
    @ingroup ISTL_Solvers
    @brief Base class for different variants of the conjugate gradient method.

    Solves a linear operator equation \f$Ax=b\f$ with \f$A:\ X\mapsto Y\f$.

    @tparam Domain type of the domain space \f$X\f$
    @tparam Range type of the range space \f$Y\f$
    @tparam BaseImpl Policy class specifiying the implementation (one of CGBase, TCGBase, RCGBase, TRCGBase in namespace CGDetail)
   */
  template <class Domain, class Range,
            template <class,class> class BaseImpl = CGDetail::CGBase>
  class CGStepImpl :
      public BaseImpl<Domain,Range> ,
      public Mixin::IterativeRefinements
  {
    using Base = BaseImpl<Domain,Range>;
    using Base::regularize;
    using Base::reset;
    using Base::adjustRegularizedResidual;
  public:
    //! type of the domain space
    using domain_type = Domain;
    //! type of the range space
    using range_type = Range;
    //! underlying field type
    using field_type = typename Domain::field_type;
    //! corresponding real type (same as real type for real spaces, differs for complex spaces)
    using real_type = typename FieldTraits<field_type>::real_type;

    /*!
      @brief Constructor.
      @param A linear operator
      @param P preconditioner
      @param sp scalar product
     */
    CGStepImpl(LinearOperator<domain_type,range_type>& A,
               Preconditioner<domain_type,range_type>& P,
               ScalarProduct<domain_type>& sp) :
      Base(sp),
      A_(A), P_(P), sp_(sp)
    {}

    /*!
      @brief Constructor using scalar product generated by SeqScalarProduct<domain_type>.
      @param A linear operator
      @param P preconditioner
     */
    CGStepImpl(LinearOperator<domain_type,range_type>& A,
               Preconditioner<domain_type,range_type>& P) :
      Base(),
      A_(A), P_(P), ssp_(), sp_(ssp_)
    {}

    /*!
      @brief Initialize conjugate gradient step.

      This method must be called before starting the conjugate gradient iteration.
      It is responsible for initializing the relevant quantities.

      @param x initial iterate
      @param b initial right hand side
     */
    void init(domain_type& x, range_type& b)
    {
      reset();
      P_.pre(x,b);

      residual_ = std::make_unique<range_type>(b);
      A_.applyscaleadd(-1,x,*residual_);
      Qr_ = std::make_unique<domain_type>(x);
      dx_ = nullptr;
      Base::init(*residual_);
    }

    /*!
      @brief Perform one step of the conjugate gradient method.
      @param x current iterate
      @param b current right hand side
     */
    void compute(domain_type& x, range_type& b)
    {
      applyPreconditioner(*residual_);
      computeSearchDirection();
      computeResidualNormWithRespectToPreconditioner();

      range_type Adx(b);
      A_.apply(*dx_,Adx);

      computeInducedStepLength(Adx);
      convexityCheck(x);
      computeStepLengthParameter();

      updateIterate(x);
      updateResidual(Adx);
    }

    /*!
      @brief Post-process final iterate, i.e. apply @code{.cpp} P_.post(x) @endcode
      @param x final iterate
     */
    void postProcess(domain_type& x)
    {
      P_.post(x);
    }

    /**
     * @brief Access scaling for the conjugate search direction.
     */
    real_type alpha() const
    {
      return alpha_;
    }

    /**
     * @brief Access length of conjugate search direction with respect to the energy norm.
     */
    real_type length() const
    {
      return dxAdx_;
    }

    /**
     * @brief Access norm of residual with respect to the norm induced by the preconditioner.
     *
     */
    real_type preconditionedResidualNorm() const
    {
      return sigma_;
    }

    /**
     * @brief Access norm of residual with respect to the employed scalar product (in general the l2-norm).
     */
    real_type residualNorm() const
    {
      return sp_.norm(*residual_);
    }

  private:
    void convexityCheck(domain_type& x)
    {
      if (dxAdx_ <= 0 )
      {
        if( this->verbosityLevel() > 1 ) std::cout << "    " << ": non-positive curvature: " << dxAdx_ << std::endl;
        this->treatNonconvexity(dxAdx_,x,*dx_);
      }
    }

    auto computeResidualNormWithRespectToPreconditioner()
    {
      using std::abs;
      sigma_ = abs( sp_.dot(*residual_,*Qr_) );
    }

    void computeInducedStepLength(const range_type& Aq)
    {
      dxAdx_ = sp_.dot(*dx_,Aq);
      regularize(dxAdx_,*dx_);
    }

    void computeStepLengthParameter()
    {
      alpha_ = sigma_/dxAdx_;
    }

    void computeSearchDirection()
    {
      if( dx_ == nullptr){
        dx_ = std::make_unique<domain_type>(*Qr_);
        return;
      }
      using std::abs;
      auto beta = abs( sp_.dot(*residual_,*Qr_) )/sigma_;
      *dx_ *= beta; *dx_ += *Qr_;
      Base::adjustPreconditionedSearchDirection(beta,*residual_);
    }

    void updateIterate(domain_type& x)
    {
      x.axpy(alpha_,*dx_);
    }

    void updateResidual(const range_type& Aq)
    {
      residual_->axpy(-alpha_,Aq);
      adjustRegularizedResidual(alpha_,*residual_);
    }

    void applyPreconditioner(range_type r) const
    {
      P_.apply(*Qr_,r);
      for(auto i=0u; i<iterativeRefinements(); ++i)
      {
        A_.applyscaleadd(-1.,*Qr_,r);
        auto dQr = range_type(r);
        P_.apply(dQr,r);
        *Qr_ += dQr;
      }
    }

    LinearOperator<domain_type,range_type>& A_;
    Preconditioner<domain_type,range_type>& P_;
    SeqScalarProduct<domain_type> ssp_;
    ScalarProduct<domain_type>& sp_;
    std::unique_ptr<range_type> residual_ = nullptr;
    std::unique_ptr<domain_type> Qr_ = nullptr, dx_ = nullptr;
    real_type alpha_ = -1, sigma_ = -1, dxAdx_ = -1;
  };
}

#endif // DUNE_CONJUGATE_GRADIENT_STEP_HH
