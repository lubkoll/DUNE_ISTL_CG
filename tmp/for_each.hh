#ifndef DUNE_TMP_FOR_EACH_HH
#define DUNE_TMP_FOR_EACH_HH

#include "dune/common/typetraits.hh"
#include "logic.hh"

namespace Dune
{
  namespace TMP
  {
    template <class,class>
    struct DefaultCombiner
    {
      template <class... Args>
      struct apply
      {
        using type = Empty;
      };
    };

    //! @cond
    namespace ForEachDetail
    {
      template <class... Args>
      struct ComposeImpl
      {
        template <class First, class Second,
                  bool = !std::is_same< typename First::template apply<Args...>::type,Empty>::value,//HasNestedType::type< typename First::template apply<Args...> > ,
                  bool = !std::is_same< typename Second::template apply<Args...>::type,Empty>::value>//HasNestedType::type< typename Second::template apply<Args...> > >
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
    }
    //! @endcond


    template <class First, class Second>
    struct Compose
    {
      template <class... Args>
      struct apply
      {
        using type = typename ForEachDetail::ComposeImpl<Args...>::template apply<First,Second>::type;
      };
    };


    template <class Operation, template <class,class> class Combiner, class... Sequence>
    struct VariadicForEach;

    template <class Operation, template <class,class> class Combiner, class Element, class... Sequence>
    struct VariadicForEach<Operation,Combiner,Element,Sequence...>
    {
      template <class... Args>
      struct apply
      {
        using type =
        typename Combiner<
          typename Operation::template apply<Element>::type,
          VariadicForEach<Operation,Combiner,Sequence...>
        >::template apply<Args...>::type;
      };
    };

    template <class Operation, template <class,class> class Combiner, class LastElement>
    struct VariadicForEach<Operation,Combiner,LastElement>
    {
      template <class... Args>
      struct apply
      {
        using type = typename Operation::template apply<LastElement>::type::template apply<Args...>::type;
      };
    };


    template <class Operation, template <class,class> class Combiner = DefaultCombiner>
    struct VariadicApply
    {
      template <class... OptionalBases>
      struct apply;

      template<class OptionalBasis, class... OtherOptionalBases>
      struct apply<OptionalBasis,OtherOptionalBases...>
      {
        using type =
        typename Combiner<
          Nullary< Operation , OptionalBasis > ,
          Nullary< VariadicApply<Operation,Combiner>, OtherOptionalBases... >
        >::template apply<>::type;
      };

      template <class OptionalBasis>
      struct apply<OptionalBasis>
      {
        using type = typename Operation::template apply<OptionalBasis>::type;
      };
    };
  }
}

#endif // DUNE_TMP_FOR_EACH_HH
