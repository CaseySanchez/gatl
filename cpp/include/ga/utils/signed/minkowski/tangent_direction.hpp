#ifndef __GA_UTILS_SIGNED_MINKOWSKI_TANGENT_DIRECTION_HPP__
#define __GA_UTILS_SIGNED_MINKOWSKI_TANGENT_DIRECTION_HPP__

namespace ga {

	// Returns the direction parameter of a given tangent.
	template<class CoefficientType, class Expression, ndims_t N>
	constexpr decltype(auto) tangent_direction(clifford_expression<CoefficientType, Expression> const &tangent, minkowski_metric_space<N> const &mtr) {
		auto const lazy = make_lazy_context(tangent);
		return lazy.eval(op(lcont(-(e(c<N + 1>) + e(c<N + 2>)), lazy.argument<0>(), mtr), e(c<N + 1>) + e(c<N + 2>), mtr));
	}

}

#endif // __GA_UTILS_SIGNED_MINKOWSKI_TANGENT_DIRECTION_HPP__