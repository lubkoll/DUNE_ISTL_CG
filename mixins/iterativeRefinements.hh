#ifndef DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH
#define DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for iterative refinements.
     */
    class IterativeRefinements : public MixinConnection<IterativeRefinements>
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
        notify();
      }

      /**
       * @brief Access number of iterative refinements.
       * @return number of iterative refinements
       */
      unsigned iterativeRefinements() const noexcept
      {
        return iterativeRefinements_;
      }

      //! Update function for a simplified observer pattern
      void update(IterativeRefinements* changed)
      {
        setIterativeRefinements( changed->iterativeRefinements() );
      }

    private:
      unsigned iterativeRefinements_;
    };
  }
}
#endif // DUNE_MIXIN_ITERATIVE_REFINEMENTS_HH
