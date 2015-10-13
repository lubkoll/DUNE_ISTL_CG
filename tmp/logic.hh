#ifndef DUNE_TMP_LOGIC_HH
#define DUNE_TMP_LOGIC_HH

#include <type_traits>

#include "empty.hh"

namespace Dune
{
  namespace TMP
  {
    struct True
    {
      template <class>
      struct apply
      {
        static constexpr bool value = true;
      };
    };

    struct False
    {
      template <class>
      struct apply
      {
        static constexpr bool value = false;
      };
    };


    template <class Type0, class Type1>
    struct And
    {
      template <class Arg>
      struct apply
      {
        static constexpr bool value = Type0::template apply<Arg>::value && Type1::template apply<Arg>::value;
      };
    };


    template <class Type0, class Type1>
    struct Or
    {
      template <class Arg>
      struct apply
      {
        static constexpr bool value = Type0::template apply<Arg>::value || Type1::template apply<Arg>::value;
      };
    };


    template <class Derived>
    struct IsDerivedFrom
    {
      template <class Base>
      struct apply
      {
        static constexpr bool value = std::is_base_of<Base,Derived>::value;
      };
    };


    template <class Derived>
    struct IsNotDerivedFrom
    {
      template <class Base>
      struct apply
      {
        static constexpr bool value = !std::is_base_of<Base,Derived>::value;
      };
    };


    template <template <class> class Unary, class Sequence>
    struct OrUnaryToSequence
    {
      template <class Arg>
      struct apply
      {
        static constexpr bool value = Unary<typename Sequence::type>::template apply<Arg>::value || OrUnaryToSequence<Unary,typename Sequence::next>::template apply<Arg>::value;
      };
    };

    template <template <class> class Unary>
    struct OrUnaryToSequence<Unary,Empty>
    {
      template <class Arg>
      struct apply
      {
        static constexpr bool value = false;
      };
    };
  }
}
#endif // DUNE_TMP_LOGIC_HH
