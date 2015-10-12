#ifndef DUNE_MIXIN_VERBOSITY_HH
#define DUNE_MIXIN_VERBOSITY_HH

#include "mixinConnection.hh"

namespace Dune
{
  namespace Mixin
  {
    /**
     * @ingroup MixinGroup
     * @brief %Mixin class for verbosity.
     */
    class Verbosity : public MixinConnection<Verbosity>
    {
    public:
      /**
       * @brief Constructor.
       * @param verbosityLevel verbosity level (0=silent,...)
       */
      explicit Verbosity(unsigned verbosityLevel = 0) noexcept
        : verbosityLevel_{verbosityLevel}
      {}

      /**
       * @brief Enable/disable verbosity.
       * @param verbose true: if verbosityLevel = 0, set verbosityLevel = 1; false: if set verbosityLevel = 0
       */
      void setVerbosity(bool verbose) noexcept
      {
        if( verbose ) verbosityLevel_ = 1;
        else verbosityLevel_ = 0;
        notify();
      }

      /**
       * @brief Check if verbosity is turned on.
       * @return true if verbosityLevel > 0
       */
      bool verbose() const noexcept
      {
        return verbosityLevel_ > 0;
      }

      /**
       * @brief Set verbosity level.
       * @param level verbosity level
       */
      void setVerbosityLevel(unsigned level) noexcept
      {
        verbosityLevel_ = level;
        notify();
      }

      /**
       * @brief Access verbosity level.
       * @return verbosity level
       */
      unsigned verbosityLevel() const noexcept
      {
        return verbosityLevel_;
      }

      //! Update function for a simplified observer pattern
      void update(Verbosity* changed)
      {
        setVerbosityLevel(changed->verbosityLevel());
      }

    private:
      unsigned verbosityLevel_;
    };
  }
}

#endif // DUNE_MIXIN_VERBOSITY_HH
