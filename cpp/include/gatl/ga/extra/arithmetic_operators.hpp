/**
Copyright (C) 2018 Leandro Augusto Frata Fernandes

author     : Fernandes, Leandro A. F.
e-mail     : laffernandes@ic.uff.br
home page  : http://www.ic.uff.br/~laffernandes
repository : https://github.com/laffernandes/gatl.git

This file is part of The Geometric Algebra Template Library (GATL).

GATL is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GATL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GATL. If not, see <https://www.gnu.org/licenses/>.
/**/

#ifndef __GA_EXTRA_ARITHMETIC_OPERATORS_HPP__
#define __GA_EXTRA_ARITHMETIC_OPERATORS_HPP__

namespace ga {

    template<typename LeftCoefficientType, typename LeftExpression, typename RightCoefficientType, typename RightCoefficient>
    constexpr decltype(auto) operator/(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, scalar_clifford_expression<RightCoefficientType, RightCoefficient> const &rhs) {
        return igp(lhs, rhs, detail::real_metric_space());
    }

    template<typename LeftCoefficientType, typename LeftExpression, typename RightType, typename = std::enable_if_t<!is_clifford_expression_v<RightType> > >
    constexpr decltype(auto) operator/(clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, RightType const &rhs) {
        return igp(lhs, scalar(rhs), detail::real_metric_space());
    }

    template<typename LeftType, typename RightCoefficientType, typename RightCoefficient, typename = std::enable_if_t<!is_clifford_expression_v<LeftType> > >
    constexpr decltype(auto) operator/(LeftType const &lhs, scalar_clifford_expression<RightCoefficientType, RightCoefficient> const &rhs) {
        return igp(scalar(lhs), rhs, detail::real_metric_space());
    }


}

#endif // __GA_EXTRA_ARITHMETIC_OPERATORS_HPP__