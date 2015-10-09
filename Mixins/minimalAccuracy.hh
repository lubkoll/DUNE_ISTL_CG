#ifndef DUNE_MIXIN_MINIMAL_ACCURACY_HH
#define DUNE_MIXIN_MINIMAL_ACCURACY_HH

#include <cassert>

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for minimal accuracy.
     */
    template <class real_type=double>
    class MinimalAccuracy
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
        minimalAccuracy_ = accuracy;
        assert(minimalAccuracy_ >= 0);
      }

      /**
       * @brief Access minimal accuracy.
       * @return minimal accuracy
       */
      real_type minimalAccuracy() const
      {
        return minimalAccuracy_;
      }

    private:
      real_type minimalAccuracy_;
    };
  }
}

#endif // DUNE_MIXIN_MINIMAL_ACCURACY_HH
