#ifndef DUNE_TMP_LOGIC_HH
#define DUNE_TMP_LOGIC_HH

#include <type_traits>

#include "dune/common/typetraits.hh"
#include "../voider.hh"

namespace Dune
{
  namespace TMP
  {
    struct True
    {
      template <class...>
      struct apply
      {
        using type = std::true_type;
      };
    };

    struct False
    {
      template <class...>
      struct apply
      {
        using type = std::false_type;
      };
    };


    template <class Operation, class OtherOperation>
    struct And
    {
      template <class... Args>
      struct apply
      {
        using type = std::integral_constant< bool , Operation::template apply<Args...>::type::value && OtherOperation::template apply<Args...>::type::value >;
      };
    };


    template <class Operation, class OtherOperation>
    struct Or
    {
      template <class... Args>
      struct apply
      {
        using type = std::integral_constant< bool , Operation::template apply<Args...>::type::value || OtherOperation::template apply<Args...>::type::value >;
      };
    };


    template <class Operation>
    struct Not
    {
      template <class... Args>
      struct apply
      {
        using type = std::integral_constant< bool , !Operation::template apply<Args...>::type::value>;
      };
    };


    template <class Derived>
    struct BaseOf
    {
      template <class Base>
      struct apply
      {
        using type = std::is_base_of<Base,Derived>;
      };
    };


    template <class Operation, class... Args>
    struct Nullary
    {
      template <class...>
      struct apply
      {
        using type = typename Operation::template apply<Args...>::type;
      };
    };


    template <class Derived>
    using NotBaseOf = Not< BaseOf<Derived> >;



    template <template <class...> class Operation>
    struct Bind
    {
      template <class... Args>
      struct apply
      {
        using type = Operation<Args...>;
      };
    };


    template <class Operation>
    struct StoreIf
    {
      template <class Type>
      struct apply
      {
        using type = typename std::conditional< Operation::template apply<Type>::type::value , Type , Empty>::type;
      };
    };


    template <class Derived>
    using StoreIfDerivedFrom = StoreIf< BaseOf<Derived> >;


    template <class Derived>
    using StoreIfNotDerivedFrom = StoreIf< NotBaseOf<Derived> >;
  }
}
#endif // DUNE_TMP_LOGIC_HH
