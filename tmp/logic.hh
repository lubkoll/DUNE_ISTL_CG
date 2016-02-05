#ifndef DUNE_TMP_LOGIC_HH
#define DUNE_TMP_LOGIC_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

namespace Dune
{
  namespace TMP
  {
    //! Template meta-function that always evaluates to std::true_type.
    struct True
    {
      template <class...>
      struct apply
      {
        using type = std::true_type;
      };
    };

    //! Template meta-function that always evaluates to std::false_type.
    struct False
    {
      template <class...>
      struct apply
      {
        using type = std::false_type;
      };
    };


    //! Convenient access to Operation::apply<Args...>::type.
    template <class Operation, class... Args>
    using Apply = typename Operation::template apply<Args...>::type;


    //! @cond
    namespace Impl
    {
      template <class Operation, class OtherOperation>
      struct And
      {
        template <class... Args>
        struct apply
        {
          using type = std::integral_constant< bool , Apply<Operation,Args...>::value && Apply<OtherOperation,Args...>::value >;
        };
      };

      template <class Operation, class OtherOperation>
      struct Or
      {
        template <class... Args>
        struct apply
        {
          using type = std::integral_constant< bool , Apply<Operation,Args...>::value || Apply<OtherOperation,Args...>::value >;
        };
      };

    }
    //! @endcond


    //! Logical "and" for meta-functions that return std::true_type or std::false_type.
    struct And
    {
      template <class First, class Second>
      struct apply
      {
        using type = Impl::And<First,Second>;
      };
    };


    //! Logical "or" for meta-functions that return std::true_type or std::false_type.
    struct Or
    {
      template <class Operation, class OtherOperation>
      struct apply
      {
        using type = Impl::Or<Operation,OtherOperation>;
      };
    };


    //! Logical "not" for meta-functions that return std::true_type or std::false_type.
    template <class Operation>
    struct Not
    {
      template <class... Args>
      struct apply
      {
        using type = std::integral_constant< bool , !Apply<Operation,Args...>::value>;
      };
    };


    //! Meta-function that checks if its argument is a base class of Derived.
    template <class Derived>
    struct BaseOf
    {
      template <class Base>
      struct apply
      {
        using type = std::is_base_of<Base,Derived>;
      };
    };


    //! Meta-function that checks if its argument is not a base class of Derived.
    template <class Derived>
    using NotBaseOf = Not< BaseOf<Derived> >;


    //! Bind Args... to Operation::apply. Thus apply becomes a nullary meta-function.
    template <class Operation, class... Args>
    struct Bind
    {
      template <class...>
      struct apply
      {
        using type = Apply<Operation,Args...>;
      };
    };


    //! Create a meta-function that takes its arguments to generate Operation<Args...>.
    template <template <class...> class Operation>
    struct BindOperation
    {
      template <class... Args>
      struct apply
      {
        using type = Operation<Args...>;
      };
    };


    //! Stores its argument Type if Operation::template apply<Type>::type::value evaluates to true, else stores a Empty.
    template <class Operation>
    struct StoreIf
    {
      template <class Type>
      struct apply
      {
        using type = typename std::conditional< Apply<Operation,Type>::value , Type , Empty>::type;
      };
    };


    template <class Derived>
    using StoreIfDerivedFrom = StoreIf< BaseOf<Derived> >;


    template <class Derived>
    using StoreIfNotDerivedFrom = StoreIf< NotBaseOf<Derived> >;
  }
}
#endif // DUNE_TMP_LOGIC_HH
