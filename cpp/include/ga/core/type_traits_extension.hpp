#ifndef __GA_CORE_TYPE_TRAITS_EXTENSION_HPP__
#define __GA_CORE_TYPE_TRAITS_EXTENSION_HPP__

namespace ga {

	namespace detail {

		// Returns true if T and any element in Rest has the same type with the same const-volatile qualifications or false otherwise.
		template<class T, class... Rest>
		constexpr bool is_any_v = std::disjunction_v<std::bool_constant<std::is_same_v<T, Rest> >...>;

		// A set of indices.
		template<std::size_t... Indices>
		struct indices {
			using next = indices<Indices..., sizeof...(Indices)>;
		};

		// Helper structure to build a set of indices.
		template<std::size_t Size>
		struct build_indices {
			using type = typename build_indices<Size - 1>::type::next;
		};

		template<>
		struct build_indices<0> {
			using type = indices<>;
		};

		template<std::size_t Size>
		using build_indices_t = typename build_indices<Size>::type;

		// Helper function to convert a tuple into a list-initialization structure.
		template<class Tuple, std::size_t... Indices>
		constexpr decltype(auto) _to_list_initialization(Tuple &&tuple, indices<Indices...>) {
			return { std::get<Indices>(std::move(tuple))... };
		}

		template<class Tuple>
		constexpr decltype(auto) to_list_initialization(Tuple &&tuple) {
			return _to_list_initialization(std::move(tuple), build_indices_t<std::tuple_size_v<std::remove_cv_t<std::remove_reference_t<Tuple> > > >());
		}

	}

}

#endif // __GA_CORE_TYPE_TRAITS_EXTENSION_HPP__
