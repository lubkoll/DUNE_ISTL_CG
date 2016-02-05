#ifndef DUNE_MIXIN_MAXSTEPS_HH
#define DUNE_MIXIN_MAXSTEPS_HH

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for maximal number of steps/iterations.
     */
    class MaxSteps : public MixinConnection<MaxSteps>
    {
    public:
      /**
       * @brief Constructor.
       * @param maxSteps maximal number of steps/iterations
       */
      explicit MaxSteps(unsigned maxSteps = 100) noexcept
        : maxSteps_{maxSteps}
      {}

      /**
       * @brief Set maximal number of steps/iterations for iterative solvers.
       * @param maxSteps maximal number of steps/iterations
       */
      void setMaxSteps(unsigned maxSteps) noexcept
      {
        maxSteps_ = maxSteps;
        notify();
      }

      /**
       * @brief Get maximal number of steps/iterations for iterative solvers.
       * @return maximal number of steps/iterations
       */
      unsigned maxSteps() const noexcept
      {
        return maxSteps_;
      }

      //! Update function for a simplified observer pattern
      void update(MaxSteps* changed)
      {
        setMaxSteps( changed->maxSteps() );
      }

    private:
      unsigned maxSteps_;
    };
  }
}

#endif // DUNE_MIXIN_MAXSTEPS_HH
