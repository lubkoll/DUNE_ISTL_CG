#ifndef DUNE_UTIL_VOIDER_HH
#define DUNE_UTIL_VOIDER_HH

namespace Dune
{
  /**
   * \cond DOCUMENT_IMPLEMENTATION_DETAILS
   */
  namespace GenericIterativeMethodDetail
  {
    /// helper to make gcc happy
    template <class...> struct voider { using type = void; };
  }
  /**
   * \endcond
   */

  /// Most fascinating type ever. Is void for all input types.
  template <class... Types>
  using void_t = typename GenericIterativeMethodDetail::voider<Types...>::type;
}

#endif // DUNE_UTIL_VOIDER_HH
