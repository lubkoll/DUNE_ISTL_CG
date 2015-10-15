#ifndef DUNE_TMP_FOR_EACH_HH
#define DUNE_TMP_FOR_EACH_HH

#include "dune/common/typetraits.hh"
#include "logic.hh"

namespace Dune
{
  namespace TMP
  {
    //! @cond
    namespace Impl
    {
      template <class... Args>
      struct Composer
      {
        template <class First, class Second,
                  bool = !std::is_same< typename First::template apply<Args...>::type,Empty>::value,
                  bool = !std::is_same< typename Second::template apply<Args...>::type,Empty>::value>
        struct apply
        {
          using type = Empty;
        };

        template <class First, class Second>
        struct apply<First,Second,true,false>
        {
          using type = typename First::template apply<Args...>::type;
        };

        template <class First, class Second>
        struct apply<First,Second,false,true>
        {
          using type = typename Second::template apply<Args...>::type;
        };

        template <class First, class Second>
        struct apply<First,Second,true,true>
        {
          struct type :
              First::template apply<Args...>::type ,
              Second::template apply<Args...>::type
          {};
        };
      };


      template <class First, class Second>
      struct Compose
      {
        template <class... Args>
        struct apply
        {
          using type = typename Composer<Args...>::template apply<First,Second>::type;
        };
      };


      template <class,class>
      struct DefaultCombiner
      {
        template <class... Args>
        struct apply
        {
          using type = Empty;
        };
      };
    }
    //! @endcond


    //! Always returns Dune::Empty.
    struct DefaultCombiner
    {
      template <class First, class Second>
      struct apply
      {
        using type = Impl::DefaultCombiner<First,Second>;
      };
    };


    /*!
      @brief Composition of return types from meta-functions.

      If both meta-functions return types different from Dune::Empty than a class that inherits from
      both results is returned.
      If both return Dune::Empty then returns Dune::Empty
      Else returns the results that differs from Dune::Empty
     */
    struct Compose
    {
      template <class First, class Second>
      struct apply
      {
        using type = Impl::Compose<First,Second>;
      };
    };


    /*!
      @tparam Operation operation to apply on each element of Sequence
      @tparam Combiner combines the individual results
     */
    template <class Operation, class Combiner, class... Sequence>
    struct VariadicForEach;

    template <class Operation, class Combiner, class Element, class... Sequence>
    struct VariadicForEach<Operation,Combiner,Element,Sequence...>
    {
      template <class... Args>
      struct apply
      {
        using type =
        typename Combiner::template apply<
          typename Operation::template apply<Element>::type,
          VariadicForEach<Operation,Combiner,Sequence...>
        >::type::template apply<Args...>::type;
      };
    };

    template <class Operation, class Combiner, class LastElement>
    struct VariadicForEach<Operation,Combiner,LastElement>
    {
      template <class... Args>
      struct apply
      {
        using type = typename Operation::template apply<LastElement>::type::template apply<Args...>::type;
      };
    };


    /*!
      @tparam Operation operation to apply on each argument to apply
      @tparam Combiner combines the individual results
     */
    template <class Operation, class Combiner = DefaultCombiner>
    struct VariadicApply
    {
      template <class... Args>
      struct apply;

      template<class Arg, class... Args>
      struct apply<Arg,Args...>
      {
        using type =
        typename Combiner::template apply<
          Bind< Operation , Arg > ,
          Bind< VariadicApply<Operation,Combiner>, Args... >
        >::type::template apply<>::type;
      };

      template <class Arg>
      struct apply<Arg> : Apply<Operation,Arg>
      {};
    };
  }
}

#endif // DUNE_TMP_FOR_EACH_HH
