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

#ifndef __GA_MODEL_SIGNED_CSTA_SPACETIME_VECTOR_HPP__
#define __GA_MODEL_SIGNED_CSTA_SPACETIME_VECTOR_HPP__

namespace ga {

    // Initializes a multivector representation of a Spacetime (STA) vector
    // using coordinates (t, x, y, z) in the base Minkowski space.
    // The result has zero conformal (ep, em) components.
    template<typename T, typename X, typename Y, typename Z>
    GA_HOST_DEVICE constexpr decltype(auto) spacetime_vector(csta_metric_space const &mtr, T &&t, X &&x, Y &&y, Z &&z) GA_NOEXCEPT {
        return vector(mtr, std::move(t), c<0>, std::move(x), std::move(y), std::move(z), c<0>);
    }

    namespace detail {

        // Helper function to adapt the iterator-based spacetime_vector().
        template<typename IteratorType, std::size_t... Indices>
        GA_ALWAYS_INLINE GA_HOST_DEVICE constexpr decltype(auto) make_spacetime_vector_using_iterator(csta_metric_space const &mtr, IteratorType begin, std::index_sequence<Indices...>) GA_NOEXCEPT {
            return spacetime_vector(mtr, *(begin + Indices)...);
        }

    }

    // Initializes a multivector representation of a Spacetime vector using
    // the given iterator to provide the four coordinates (t, x, y, z).
    template<typename IteratorType, std::enable_if_t<detail::is_iterator_v<IteratorType>, int> = 0>
    GA_HOST_DEVICE constexpr decltype(auto) spacetime_vector(csta_metric_space const &mtr, IteratorType begin, IteratorType end) GA_NOEXCEPT {
        assert(4 == std::distance(begin, end));
        return detail::make_spacetime_vector_using_iterator(mtr, begin, std::make_index_sequence<4>{});
    }

}

#endif // __GA_MODEL_SIGNED_CSTA_SPACETIME_VECTOR_HPP__
