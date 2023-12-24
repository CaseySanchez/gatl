/* Copyright (C) Leandro Augusto Frata Fernandes
 * 
 * author     : Fernandes, Leandro A. F.
 * e-mail     : laffernandes@ic.uff.br
 * home page  : http://www.ic.uff.br/~laffernandes
 * repository : https://github.com/laffernandes/gatl.git
 * 
 * This file is part of The Geometric Algebra Template Library (GATL).
 * 
 * GATL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * GATL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GATL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef __GA_EXTRA_INVERSE_GEOMETRIC_PRODUCT_HPP__
#define __GA_EXTRA_INVERSE_GEOMETRIC_PRODUCT_HPP__

namespace ga {

    template<typename LeftCoefficientType, typename LeftExpression, typename RightCoefficientType, typename RightExpression, typename MetricSpaceType>
    constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs, metric_space<MetricSpaceType> const &mtr) {
        auto [lazy, lhs_, rhs_] = make_lazy_context_tuple(lhs, rhs);
        return lazy.eval(gp(lhs_, inv(rhs_, mtr), mtr));
    }

    template<typename LeftCoefficientType, typename LeftExpression, typename RightCoefficientType, typename RightCoefficient>
    constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, scalar_clifford_expression<RightCoefficientType, RightCoefficient> const &rhs) {
        return igp(lhs, rhs, detail::real_metric_space());
    }

    template<typename LeftCoefficientType, typename LeftExpression, typename MetricSpaceType>
    constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, NonCliffordExpressionType auto const &rhs, metric_space<MetricSpaceType> const &) {
        return igp(lhs, scalar(rhs), detail::real_metric_space());
    }

    template<typename LeftCoefficientType, typename LeftExpression>
    constexpr decltype(auto) igp(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, NonCliffordExpressionType auto const &rhs) {
        return igp(lhs, scalar(rhs), detail::real_metric_space());
    }

    template<typename RightCoefficientType, typename RightExpression, typename MetricSpaceType>
    constexpr decltype(auto) igp(NonCliffordExpressionType auto const &lhs, clifford_expression<RightCoefficientType, RightExpression> const &rhs, metric_space<MetricSpaceType> const &mtr) {
        return igp(scalar(lhs), rhs, mtr);
    }

    template<typename RightCoefficientType, typename RightCoefficient>
    constexpr decltype(auto) igp(NonCliffordExpressionType auto const &lhs, scalar_clifford_expression<RightCoefficientType, RightCoefficient> const &rhs) {
        return igp(scalar(lhs), rhs, detail::real_metric_space());
    }

    template<typename MetricSpaceType>
    constexpr decltype(auto) igp(NonCliffordExpressionType auto const &lhs, NonCliffordExpressionType auto const &rhs, metric_space<MetricSpaceType> const &) {
        return igp(scalar(lhs), scalar(rhs), detail::real_metric_space());
    }

    constexpr decltype(auto) igp(NonCliffordExpressionType auto const &lhs, NonCliffordExpressionType auto const &rhs) {
        return igp(scalar(lhs), scalar(rhs), detail::real_metric_space());
    }

}

#endif // __GA_EXTRA_INVERSE_GEOMETRIC_PRODUCT_HPP__
