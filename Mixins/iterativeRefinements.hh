#ifndef DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH
#define DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for iterative refinements.
     */
    class IterativeRefinements
    {
    public:
      /**
       * @brief Constructor.
       * @param refinements number of iterative refinements.
       */
      explicit IterativeRefinements(unsigned refinements = 0) noexcept
        : iterativeRefinements_{refinements}
      {}

      /**
       * @brief Set number of iterative refinements.
       * @param refinements number of iterative refinements
       */
      void setIterativeRefinements(unsigned refinements) noexcept
      {
        iterativeRefinements_ = refinements;
      }

      /**
       * @brief Access number of iterative refinements.
       * @return number of iterative refinements
       */
      unsigned iterativeRefinements() const noexcept
      {
        return iterativeRefinements_;
      }

    private:
      unsigned iterativeRefinements_;
    };
  }
}
#endif // DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH
