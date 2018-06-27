#ifndef __GA2H_HPP__
#define __GA2H_HPP__

#include <ga/core.hpp>
#include <ga/extra.hpp>
#include <ga/utils/homogeneous.hpp>

namespace ga2h {

	using namespace ga;

	_GA_UTILS_HOMOGENEOUS_ALGEBRA_DEFINITION(space, basis_vectors_names, 2, "e1", "e2")

	static auto const e1 = e(c<1>);
	static auto const e2 = e(c<2>);

	_GA_CORE_OVERLOAD(space)
	_GA_EXTRA_OVERLOAD(space, basis_vectors_names)

}

#endif // __GA2H_HPP__
