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

#ifndef __GA_MODEL_SIGNED_CSTA_MACRO_FOR_ALGEBRA_DEFINITION_HPP__
#define __GA_MODEL_SIGNED_CSTA_MACRO_FOR_ALGEBRA_DEFINITION_HPP__

// Macro to define the Conformal Spacetime Algebra (CSTA), Cl(2,4).
//
// Basis vectors (in signed_metric_space<2,4,0> index order):
//   e(c<1>) = et  (+1, timelike)
//   e(c<2>) = ep  (+1, conformal positive)
//   e(c<3>) = ex  (-1, spacelike x)
//   e(c<4>) = ey  (-1, spacelike y)
//   e(c<5>) = ez  (-1, spacelike z)
//   e(c<6>) = em  (-1, conformal negative)
//
// Null conformal vectors:
//   no = (em - ep) / 2   (conformal origin)
//   ni = ep + em         (conformal infinity)
//
// STA sub-pseudoscalar:
//   Ista = rcont(I, op(ep, em, SPACE), SPACE)
//        = et ^ ex ^ ey ^ ez  (pseudoscalar of the Minkowski base space)
#define _GA_CSTA_ALGEBRA_DEFINITION(SPACE, BASIS_VECTORS_NAMES) \
    using space_t = csta_metric_space; \
    \
    static space_t const SPACE; \
    static std::array<std::string, 6> const BASIS_VECTORS_NAMES = { "et", "ep", "ex", "ey", "ez", "em" }; \
    \
    static auto const et = e(c<1>); \
    static auto const ep = e(c<2>); \
    static auto const ex = e(c<3>); \
    static auto const ey = e(c<4>); \
    static auto const ez = e(c<5>); \
    static auto const em = e(c<6>); \
    \
    static auto const no = (em - ep) / c<2>; \
    static auto const ni = ep + em; \
    \
    static auto const _0 = c<0>; \
    static auto const _1 = c<1>; \
    static auto const _2 = c<2>; \
    \
    static auto const I = pseudoscalar(SPACE); \
    static auto const Ista = rcont(I, op(ep, em, SPACE), SPACE);

#endif // __GA_MODEL_SIGNED_CSTA_MACRO_FOR_ALGEBRA_DEFINITION_HPP__
