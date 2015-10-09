#ifndef DUNE_MIXIN_MAXSTEPS_HH
#define DUNE_MIXIN_MAXSTEPS_HH

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for maximal number of steps/iterations.
     */
    class MaxSteps
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
      }

      /**
       * @brief Get maximal number of steps/iterations for iterative solvers.
       * @return maximal number of steps/iterations
       */
      unsigned maxSteps() const noexcept
      {
        return maxSteps_;
      }

    private:
      unsigned maxSteps_;
    };
  }
}

#endif // DUNE_MIXIN_MAXSTEPS_HH