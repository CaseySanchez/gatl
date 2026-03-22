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

#ifndef __GA_MODEL_SIGNED_CSTA_METRIC_SPACE_HPP__
#define __GA_MODEL_SIGNED_CSTA_METRIC_SPACE_HPP__

namespace ga {

    // Conformal Spacetime Algebra (CSTA) metric space, Cl(2,4).
    //
    // The CSTA is the conformal extension of the Spacetime Algebra (STA), Cl(1,3).
    // Adding one positive (ep) and one negative (em) conformal dimension to STA
    // yields Cl(2,4) with basis ordering:
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
    class csta_metric_space : public signed_metric_space<2, 4, 0> {
    };

}

#endif // __GA_MODEL_SIGNED_CSTA_METRIC_SPACE_HPP__
