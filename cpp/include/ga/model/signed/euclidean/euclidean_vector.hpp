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

#ifndef __GA_MODEL_SIGNED_EUCLIDEAN_EUCLIDEAN_VECTOR_HPP__
#define __GA_MODEL_SIGNED_EUCLIDEAN_EUCLIDEAN_VECTOR_HPP__

namespace ga {

	// Initializes a multivector representation of a Euclidean vector using the given coordinates expressed in the base space.
	template<ndims_t N, typename... Types>
	constexpr decltype(auto) euclidean_vector(euclidean_metric_space<N> const &mtr, Types &&... coords) {
		return vector(mtr, std::move(coords)...);
	}

}

#endif // __GA_MODEL_SIGNED_EUCLIDEAN_EUCLIDEAN_VECTOR_HPP__
