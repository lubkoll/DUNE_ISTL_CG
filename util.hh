#ifndef DUNE_UTIL_HH
#define DUNE_UTIL_HH

#include "voider.hh"

namespace Dune
{
  //! @cond
  namespace UtilDetail{ template <class,class = void> class Real_t; }
  //! @endcond


  //! Convenient access to Type::field_type.
  template <class Type>
  using field_t = typename Type::field_type;


  //! Convenient access to Type::real_type, resp. FieldTraits<Type::field_type>::real_type.
  template <class Type>
  using real_t = typename UtilDetail::Real_t<Type>::type;


  //! Convenient access to Type::type;
  template <class Type>
  using access_t = typename Type::type;


  //! @cond
  namespace UtilDetail
  {
    template <class Type, class>
    struct Real_t
    {
      using type = typename FieldTraits<Type>::real_type;
    };

    template <class Type>
    struct Real_t< Type , void_t< field_t<Type> > >
    {
      using type = typename Real_t< field_t<Type> >::type;
    };
  }
  //! @endcond
}

#endif // DUNE_UTIL_HH
