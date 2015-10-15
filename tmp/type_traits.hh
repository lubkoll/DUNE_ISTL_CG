#ifndef DUNE_TMP_TYPE_TRAITS_HH
#define DUNE_TMP_TYPE_TRAITS_HH

#include <type_traits>

#include "../voider.hh"

namespace Dune
{
  namespace TryNestedType
  {
    template <class Type>
    using type = typename Type::type;
  }

  namespace HasNestedType
  {
    template <class Type, class = void>
    struct type : std::false_type {};

    template <class Type>
    struct type< Type , void_t< TryNestedType::type<Type> > > : std::true_type {};
  }

  namespace TMP
  {
    template <class Type>
    using enable_if_t = void_t< TryNestedType::type<Type> >;
  }
}

#endif // DUNE_TMP_TYPE_TRAITS_HH
