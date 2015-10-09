#ifndef DUNE_MIXIN_RELATIVE_ACCURACY_HH
#define DUNE_MIXIN_RELATIVE_ACCURACY_HH

#include <cassert>
#include <limits>

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for relative accuracy.
     */
    template <class real_type=double>
    class RelativeAccuracy
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
        relativeAccuracy_ = accuracy;
        assert(relativeAccuracy_ >= 0);
      }

      /**
       * @brief Access relative accuracy.
       * @return relative accuracy
       */
      real_type relativeAccuracy() const
      {
        return relativeAccuracy_;
      }

    private:
      real_type relativeAccuracy_;
    };
  }
}

#endif // DUNE_MIXIN_RELATIVE_ACCURACY_HH
