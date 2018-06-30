#ifndef __GA_UTILS_SIGNED_HOMOGENEOUS_FLAT_UNIT_SUPPORT_POINT_HPP__
#define __GA_UTILS_SIGNED_HOMOGENEOUS_FLAT_UNIT_SUPPORT_POINT_HPP__

namespace ga {

	// Returns the unit support point parameter of a given k-flat.
	template<class CoefficientType, class Expression, ndims_t N>
	constexpr decltype(auto) flat_unit_support_point(clifford_expression<CoefficientType, Expression> const &flat, homogeneous_metric_space<N> const &mtr) {
		auto const lazy = make_lazy_context(flat);
		return lazy.eval(gp(lazy.argument<0>(), inv(lcont(e(c<N + 1>), lazy.argument<0>(), mtr), mtr), mtr));
	}

}

#endif // __GA_UTILS_SIGNED_HOMOGENEOUS_FLAT_UNIT_SUPPORT_POINT_HPP__