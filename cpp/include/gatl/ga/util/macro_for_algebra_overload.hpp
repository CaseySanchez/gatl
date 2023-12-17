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

#ifndef __GA_UTIL_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__
#define __GA_UTIL_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__

#define _GA_UTIL_OVERLOAD(SPACE) \
    constexpr decltype(auto) pseudoscalar() GA_NOEXCEPT { \
        return pseudoscalar(SPACE); \
    } \
    \
    template<typename... Types> requires (std::disjunction_v<std::bool_constant<!detail::is_iterator_v<Types> >...>) \
    constexpr decltype(auto) vector(Types &&... coords) GA_NOEXCEPT { \
        return ga::vector(SPACE, std::move(coords)...); \
    } \
    \
    template<typename IteratorType> requires (detail::is_iterator_v<IteratorType>) \
    constexpr decltype(auto) vector(IteratorType begin, IteratorType end) GA_NOEXCEPT { \
        return ga::vector(SPACE, begin, end); \
    }

#endif // __GA_UTIL_MACRO_FOR_ALGEBRA_OVERLOAD_HPP__
