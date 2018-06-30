#ifndef __GA2C_HPP__
#define __GA2C_HPP__

#include <ga/core.hpp>
#include <ga/extra.hpp>
#include <ga/utils/conformal.hpp>

namespace ga2c {

	using namespace ga;

	_GA_UTILS_CONFORMAL_ALGEBRA_DEFINITION(space, basis_vectors_names, 2, "e1", "e2")
		
	static auto const e1 = e(c<1>);
	static auto const e2 = e(c<2>);

	_GA_CORE_OVERLOAD(space)
	_GA_EXTRA_OVERLOAD(space, basis_vectors_names)
	_GA_UTILS_CONFORMAL_ALGEBRA_OVERLOAD(space)

}

#endif // __GA2C_HPP__
