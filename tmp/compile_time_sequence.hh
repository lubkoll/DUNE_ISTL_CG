#ifndef DUNE_TMP_COMPILE_TIME_SEQUENCE_HH
#define DUNE_TMP_COMPILE_TIME_SEQUENCE_HH

#include "empty.hh"

namespace Dune
{
  namespace TMP
  {
    template <class...>
    struct Sequence;

    template <class Arg, class... Args>
    struct Sequence<Arg,Args...>
    {
      using type = Arg;
      using next = Sequence<Args...>;
    };

    template <class Arg>
    struct Sequence<Arg>
    {
      using type = Arg;
      using next = Empty;
    };
  }
}

#endif // DUNE_TMP_COMPILE_TIME_SEQUENCE_HH
