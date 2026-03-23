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

#ifndef __GA_MODEL_SIGNED_CSTA_POINT_HPP__
#define __GA_MODEL_SIGNED_CSTA_POINT_HPP__

namespace ga {

    // Initializes a multivector representation of a conformal spacetime point
    // from coordinates (t, x, y, z) in Minkowski spacetime.
    //
    // The conformal embedding maps X = t*et + x*ex + y*ey + z*ez to:
    //   P = X + ((X^2 - 1)/2)*ep + ((X^2 + 1)/2)*em
    // where X^2 = t^2 - x^2 - y^2 - z^2 is the Minkowski inner product.
    //
    // Basis order in vector(mtr, ...):
    //   (et_coeff, ep_coeff, ex_coeff, ey_coeff, ez_coeff, em_coeff)
    template<typename T, typename X, typename Y, typename Z>
    GA_HOST_DEVICE constexpr decltype(auto) point(csta_metric_space const &mtr, T &&t, X &&x, Y &&y, Z &&z) GA_NOEXCEPT {
        auto aux = t * t - x * x - y * y - z * z;
        return vector(mtr, std::move(t), (aux - c<1>) / c<2>, std::move(x), std::move(y), std::move(z), (aux + c<1>) / c<2>);
    }

    namespace detail {

        // Helper function to adapt the iterator-based point().
        template<typename IteratorType, std::size_t... Indices>
        GA_ALWAYS_INLINE GA_HOST_DEVICE constexpr decltype(auto) make_csta_point_using_iterator(csta_metric_space const &mtr, IteratorType begin, std::index_sequence<Indices...>) GA_NOEXCEPT {
            return point(mtr, *(begin + Indices)...);
        }

    }

    // Initializes a multivector representation of a conformal spacetime point
    // using the given iterator to provide the four coordinates (t, x, y, z).
    template<typename IteratorType, std::enable_if_t<detail::is_iterator_v<IteratorType>, int> = 0>
    GA_HOST_DEVICE constexpr decltype(auto) point(csta_metric_space const &mtr, IteratorType begin, IteratorType end) GA_NOEXCEPT {
        assert(4 == std::distance(begin, end));
        return detail::make_csta_point_using_iterator(mtr, begin, std::make_index_sequence<4>{});
    }

}

#endif // __GA_MODEL_SIGNED_CSTA_POINT_HPP__
