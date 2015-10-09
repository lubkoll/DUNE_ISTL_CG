#ifndef DUNE_MIXIN_EPS_HH
#define DUNE_MIXIN_EPS_HH

#include <cassert>
#include <cmath>
#include <limits>

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for maximal attainable accuracy \f$\varepsilon\f$.
     */
    template <class real_type=double>
    class Eps
    {
    public:
      /**
       * @brief Constructor.
       * @param eps maximal attainable accuracy \f$\varepsilon\f$
       */
      explicit Eps( real_type eps = std::numeric_limits<real_type>::epsilon() )
        : eps_{eps}
      {
        assert(eps_ > 0);
      }

      /**
       * @brief Set maximal attainable accuracy \f$\varepsilon\f$.
       * @param eps new maximal attainable accuracy
       */
      void setEps(real_type eps)
      {
        eps_ = eps;
        assert(eps_ > 0);
      }

      /**
       * @brief Access maximal attainable accuracy.
       * @return \f$\varepsilon\f$
       */
      real_type eps() const
      {
        return eps_;
      }

      /**
       * @brief Access square root of maximal attainable accuracy.
       * @return \f$\sqrt\varepsilon\f$
       */
      real_type sqrtEps() const
      {
        return sqrt(eps_);
      }

      /**
       * @brief Get third root of maximal attainable accuracy.
       * @return \f$\varepsilon^{1/3}\f$
       */
      real_type cbrtEps() const
      {
        return cbrt(eps_);
      }

    private:
      real_type eps_;
    };
  }
}

#endif // DUNE_MIXIN_EPS_HH
