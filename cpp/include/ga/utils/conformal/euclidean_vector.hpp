#ifndef __GA_UTILS_CONFORMAL_EUCLIDEAN_VECTOR_HPP__
#define __GA_UTILS_CONFORMAL_EUCLIDEAN_VECTOR_HPP__

namespace ga {

	// Initializes a multivector representation of a Euclidean vector using the given coordinates expressed in the base space.
	template<ndims_t N, class... Types>
	constexpr decltype(auto) euclidean_vector(conformal_metric_space<N> const &mtr, Types &&... coords) {
		return vector(mtr, std::move(coords)..., c<0>, c<0>);
	}

}

#endif // __GA_UTILS_CONFORMAL_EUCLIDEAN_VECTOR_HPP__