#ifndef DUNE_OPTIONAL_BASE_CLASS_HH
#define DUNE_OPTIONAL_BASE_CLASS_HH

#include <type_traits>

#include "empty.hh"

namespace Dune
{
  namespace TMP
  {
    template <class First, class Second,
              bool = !std::is_same<First,Empty>::value,
              bool = !std::is_same<Second,Empty>::value>
    struct ComposeClass
    {
      using type = Empty;
    };

    template <class First, class Second>
    struct ComposeClass<First,Second,true,false>
    {
      using type = First;
    };

    template <class First, class Second>
    struct ComposeClass<First,Second,false,true>
    {
      using type = Second;
    };

    template <class First, class Second>
    struct ComposeClass<First,Second,true,true>
    {
      struct type : First , Second
      {};
    };


    template <class Decider, class... OptionalBases>
    struct BaseClassesIf;

    template <class Decider, class OptionalBasis, class... OtherOptionalBases>
    struct BaseClassesIf<Decider,OptionalBasis,OtherOptionalBases...>
    {
      using type =
      typename ComposeClass<
        typename std::conditional< Decider::template apply<OptionalBasis>::value , OptionalBasis , Empty >::type,
        typename BaseClassesIf< Decider , OtherOptionalBases... >::type
      >::type;
    };

    template <class Decider, class OptionalBasis>
    struct BaseClassesIf<Decider,OptionalBasis>
    {
      using type = typename std::conditional< Decider::template apply<OptionalBasis>::value , OptionalBasis , Empty >::type;
    };
  }
}

#endif // DUNE_OPTIONAL_BASE_CLASS_HH
