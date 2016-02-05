#ifndef DUNE_MIXIN_ABSOLUTE_ACCURACY_HH
#define DUNE_MIXIN_ABSOLUTE_ACCURACY_HH

#include <cassert>
#include <limits>

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for absolute accuracy.
     */
    template <class real_type=double>
    class AbsoluteAccuracy : public MixinConnection< AbsoluteAccuracy<real_type> >
    {
    public:
      /**
       * @brief Constructor.
       * @param accuracy absolute accuracy
       */
      explicit AbsoluteAccuracy(real_type accuracy = std::numeric_limits<real_type>::epsilon())
        : absoluteAccuracy_{accuracy}
      {
        assert(absoluteAccuracy_ >= 0);
      }

      /**
       * @brief Set absolute accuracy.
       * @param accuracy absolute accuracy
       */
      void setAbsoluteAccuracy(real_type accuracy)
      {
        assert(accuracy >= 0);
        absoluteAccuracy_ = accuracy;
        this->notify();
      }

      /**
       * @brief Access absolute accuracy.
       * @return absolute accuracy
       */
      real_type absoluteAccuracy() const
      {
        return absoluteAccuracy_;
      }

      //! Update function for a simplified observer pattern with MixinConnection.
      void update(AbsoluteAccuracy* changed)
      {
        setAbsoluteAccuracy( changed->absoluteAccuracy() );
      }

    private:
      real_type absoluteAccuracy_;
    };
  }
}

#endif // DUNE_MIXIN_ABSOLUTE_ACCURACY_HH
