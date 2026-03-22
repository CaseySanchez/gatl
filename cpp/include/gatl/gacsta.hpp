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

#ifndef __GACSTA_HPP__
#define __GACSTA_HPP__

#include "ga/core.hpp"
#include "ga/util.hpp"
#include "ga/extra.hpp"
#include "ga/model/csta.hpp"

// Geometry functions require ga/extra.hpp (for operator/, igp, inv) to be
// included first. They are included here, after ga/extra.hpp, rather than
// from ga/model/csta.hpp, for the same reason the similar functions in the
// minkowski and conformal models are only available after ga/extra.hpp.
#include "ga/model/signed/csta/spacetime_vector.hpp"
#include "ga/model/signed/csta/point.hpp"
#include "ga/model/signed/csta/flat_direction.hpp"
#include "ga/model/signed/csta/flat_location.hpp"
#include "ga/model/signed/csta/round_direction.hpp"
#include "ga/model/signed/csta/round_location.hpp"
#include "ga/model/signed/csta/round_size_sqr.hpp"
#include "ga/model/signed/csta/tangent_direction.hpp"
#include "ga/model/signed/csta/tangent_location.hpp"

// Conformal Spacetime Algebra (CSTA): Cl(2,4)
//
// The CSTA is the conformal extension of the Spacetime Algebra (STA), Cl(1,3).
// It embeds Minkowski spacetime in a 6-dimensional space with signature (2,4,0).
//
// Basis vectors (east-coast / physics convention, time-positive):
//   et  — timelike basis vector,          et^2 = +1
//   ep  — conformal positive null-basis,  ep^2 = +1
//   ex  — spacelike x basis vector,       ex^2 = -1
//   ey  — spacelike y basis vector,       ey^2 = -1
//   ez  — spacelike z basis vector,       ez^2 = -1
//   em  — conformal negative null-basis,  em^2 = -1
//
// Derived null conformal vectors:
//   no = (em - ep) / 2   (conformal origin)
//   ni = ep + em         (conformal infinity / point at infinity)
//
// Pseudoscalars:
//   I    = et ^ ep ^ ex ^ ey ^ ez ^ em   (full CSTA pseudoscalar)
//   Ista = et ^ ex ^ ey ^ ez             (STA sub-pseudoscalar)
//
// Conformal point embedding for X = (t, x, y, z) in Minkowski spacetime:
//   P = t*et + x*ex + y*ey + z*ez + ((t^2-x^2-y^2-z^2-1)/2)*ep
//                                  + ((t^2-x^2-y^2-z^2+1)/2)*em
//   P = spacetime_vector(t,x,y,z) + (X^2/2)*ni + no
//   where X^2 = t^2 - x^2 - y^2 - z^2 is the Minkowski inner product.

namespace gacsta {

    using namespace ga;

    _GA_CSTA_ALGEBRA_DEFINITION(space, basis_vectors_names)

    _GA_CORE_OVERLOAD(space)
    _GA_UTIL_OVERLOAD(space)
    _GA_EXTRA_OVERLOAD(space, basis_vectors_names)
    _GA_CSTA_ALGEBRA_OVERLOAD(space)

}

#endif // __GACSTA_HPP__
