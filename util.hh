#ifndef DUNE_UTIL_HH
#define DUNE_UTIL_HH

#include "voider.hh"

namespace Dune
{
  //! @cond
  namespace UtilDetail
  {
    template <class,class = void> struct Real_t;
    class Empty {};
  }
  //! @endcond

  //! Convenient access to Type::field_type.
  template <class Type>
  using field_t = typename Type::field_type;


  //! Convenient access to Type::real_type, resp. FieldTraits<Type::field_type>::real_type.
  template <class Type>
  using real_t = typename UtilDetail::Real_t<Type>::type;


  //! Convenient access to Type::type
  template <class Type>
  using access_t = typename Type::type;


  //! Evaluates to Base if std::is_base_of<Base,Derived>::value==false, else evaluates to an empty class.
  template <class Derived, class Base>
  using DeriveFromBaseIfNotAmbiguous = access_t<std::conditional< std::is_base_of<Base,Derived>::value , UtilDetail::Empty , Base > >;


  /*! @cond */
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
  /*! @endcond */
}

#endif // DUNE_UTIL_HH
