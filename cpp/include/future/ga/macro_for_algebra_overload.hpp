#ifndef __FUTURE_GA_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__
#define __FUTURE_GA_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__

#define GA_SIGNED_ALGEBRA_OVERLOAD(NAMESPACE_MNEMONIC, P, Q) \
	namespace NAMESPACE_MNEMONIC { \
		\
		using namespace future::ga; \
		\
		static signed_metric_space<P, Q> const space; \
		\
		constexpr decltype(auto) pseudoscalar() { \
			return pseudoscalar(space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) gp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return gp(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) lcont(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return lcont(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) op(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return op(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression > \
		constexpr decltype(auto) operator^(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return op(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) rcont(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return rcont(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) scp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return scp(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return igp(lhs, rhs, space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) inv(clifford_expression<CoefficientType, Expression> const &arg) { \
			return inv(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) rnorm_sqr(clifford_expression<CoefficientType, Expression> const &arg) { \
			return rnorm_sqr(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) rnorm(clifford_expression<CoefficientType, Expression> const &arg) { \
			return rnorm(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression, class PseudoscalarType, class PseudoscalarCoefficientType, class PseudoscalarExpression> \
		constexpr decltype(auto) dual(clifford_expression<CoefficientType, Expression> const &arg, clifford_expression<PseudoscalarCoefficientType, PseudoscalarExpression> const &pseudoscalar) { \
			return dual(arg, pseudoscalar, space); \
		} \
		\
		template<class Type> \
		constexpr decltype(auto) dual(Type const &arg) { \
			return dual(arg, pseudoscalar(space), space); \
		} \
		\
		template<class CoefficientType, class Expression, class PseudoscalarType, class PseudoscalarCoefficientType, class PseudoscalarExpression> \
		constexpr decltype(auto) undual(clifford_expression<CoefficientType, Expression> const &arg, clifford_expression<PseudoscalarCoefficientType, PseudoscalarExpression> const &pseudoscalar) { \
			return undual(arg, pseudoscalar, space); \
		} \
		\
		template<class Type> \
		constexpr decltype(auto) undual(Type const &arg) { \
			return undual(arg, pseudoscalar(space), space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) exp(clifford_expression<CoefficientType, Expression> const &arg) { \
			return exp(arg, space); \
		} \
		\
	}

#define GA_CONFORMAL_ALGEBRA_OVERLOAD(NAMESPACE_MNEMONIC, N) \
	namespace NAMESPACE_MNEMONIC { \
		\
		using namespace future::ga; \
		\
		static conformal_metric_space<N> const space; \
		\
		constexpr decltype(auto) pseudoscalar() { \
			return pseudoscalar(space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) gp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return gp(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) lcont(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return lcont(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) op(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return op(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) rcont(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return rcont(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) scp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return scp(lhs, rhs, space); \
		} \
		\
		template<class LeftCoefficientType, class LeftExpression, class RightCoefficientType, class RightExpression> \
		constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs) { \
			return igp(lhs, rhs, space); \
		} \
		\
				template<class CoefficientType, class Expression> \
		constexpr decltype(auto) inv(clifford_expression<CoefficientType, Expression> const &arg) { \
			return inv(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) rnorm_sqr(clifford_expression<CoefficientType, Expression> const &arg) { \
			return rnorm_sqr(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) rnorm(clifford_expression<CoefficientType, Expression> const &arg) { \
			return rnorm(arg, space); \
		} \
		\
		template<class CoefficientType, class Expression, class PseudoscalarType, class PseudoscalarCoefficientType, class PseudoscalarExpression> \
		constexpr decltype(auto) dual(clifford_expression<CoefficientType, Expression> const &arg, clifford_expression<PseudoscalarCoefficientType, PseudoscalarExpression> const &pseudoscalar) { \
			return dual(arg, pseudoscalar, space); \
		} \
		\
		template<class Type> \
		constexpr decltype(auto) dual(Type const &arg) { \
			return dual(arg, pseudoscalar(space), space); \
		} \
		\
		template<class CoefficientType, class Expression, class PseudoscalarType, class PseudoscalarCoefficientType, class PseudoscalarExpression> \
		constexpr decltype(auto) undual(clifford_expression<CoefficientType, Expression> const &arg, clifford_expression<PseudoscalarCoefficientType, PseudoscalarExpression> const &pseudoscalar) { \
			return undual(arg, pseudoscalar, space); \
		} \
		\
		template<class Type> \
		constexpr decltype(auto) undual(Type const &arg) { \
			return undual(arg, pseudoscalar(space), space); \
		} \
		\
		template<class CoefficientType, class Expression> \
		constexpr decltype(auto) exp(clifford_expression<CoefficientType, Expression> const &arg) { \
			return exp(arg, space); \
		} \
		\
	}

#endif // __FUTURE_GA_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__