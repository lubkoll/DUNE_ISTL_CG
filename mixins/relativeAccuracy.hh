#ifndef DUNE_MIXIN_RELATIVE_ACCURACY_HH
#define DUNE_MIXIN_RELATIVE_ACCURACY_HH

#include <cassert>
#include <limits>

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for relative accuracy.
     */
    template <class real_type=double>
    class RelativeAccuracy : public MixinConnection< RelativeAccuracy<real_type> >
    {
    public:
      /**
       * @brief Constructor.
       * @param accuracy relative accuracy.
       */
      explicit RelativeAccuracy(real_type accuracy = std::numeric_limits<real_type>::epsilon())
        : relativeAccuracy_{accuracy}
      {
        assert(relativeAccuracy_ >= 0);
      }

      /**
       * @brief Set relative accuracy.
       * @param accuracy relative accuracy
       */
      void setRelativeAccuracy(real_type accuracy)
      {
        assert(accuracy>= 0);
        relativeAccuracy_ = accuracy;
        this->notify();
      }

      /**
       * @brief Access relative accuracy.
       * @return relative accuracy
       */
      real_type relativeAccuracy() const
      {
        return relativeAccuracy_;
      }

      //! Update function for a simplified observer pattern
      void update(RelativeAccuracy* changed)
      {
        setRelativeAccuracy( changed->relativeAccuracy() );
      }

    private:
      real_type relativeAccuracy_;
    };
  }
}

#endif // DUNE_MIXIN_RELATIVE_ACCURACY_HH
