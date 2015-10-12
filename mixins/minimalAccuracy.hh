#ifndef DUNE_MIXIN_MINIMAL_ACCURACY_HH
#define DUNE_MIXIN_MINIMAL_ACCURACY_HH

#include <cassert>

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for minimal accuracy.
     */
    template <class real_type=double>
    class MinimalAccuracy : public MixinConnection< MinimalAccuracy<real_type> >
    {
    public:
      /**
       * @brief Constructor.
       * @param accuracy minimal accuracy
       */
      explicit MinimalAccuracy(real_type accuracy = 0.25)
        : minimalAccuracy_{accuracy}
      {
        assert(minimalAccuracy_ >= 0);
      }

      /**
       * @brief Set minimal accuracy.
       * @param accuracy minimal accuracy
       */
      void setMinimalAccuracy(real_type accuracy)
      {
        assert(accuracy >= 0);
        minimalAccuracy_ = accuracy;
        this->notify();
      }

      /**
       * @brief Access minimal accuracy.
       * @return minimal accuracy
       */
      real_type minimalAccuracy() const
      {
        return minimalAccuracy_;
      }

      //! Update function for a simplified observer pattern
      void update(MinimalAccuracy* changed)
      {
        setMinimalAccuracy( changed->minimalAccuracy() );
      }

    private:
      real_type minimalAccuracy_;
    };
  }
}

#endif // DUNE_MIXIN_MINIMAL_ACCURACY_HH
