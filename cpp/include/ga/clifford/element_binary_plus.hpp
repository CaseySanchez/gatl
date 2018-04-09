#ifndef __GA_CLIFFORD_ELEMENT_BINARY_PLUS_HPP__
#define __GA_CLIFFORD_ELEMENT_BINARY_PLUS_HPP__

namespace ga {

	namespace clifford {

		namespace detail {

			template<class LeftCoefficientType, default_bitset_t BasisBlade, class RightCoefficientType>
			constexpr decltype(auto) binary_plus_element(component<LeftCoefficientType, cbasis_blade<BasisBlade> > const &lhs, component<RightCoefficientType, cbasis_blade<BasisBlade> > const &rhs) {
				return make_component(lhs.coefficient() + rhs.coefficient(), cbasis_blade<BasisBlade>());
			}

			template<class LeftCoefficientType, default_bitset_t PossibleGrades, class RightCoefficientType>
			constexpr decltype(auto) binary_plus_element(component<LeftCoefficientType, dbasis_blade<PossibleGrades> > const &lhs, component<RightCoefficientType, dbasis_blade<PossibleGrades> > const &rhs) {
				//TODO LAZY O de cima � o original
				/**
				typedef typename std::common_type<LeftCoefficientType, RightCoefficientType, decltype(LeftCoefficientType() + RightCoefficientType())>::type coefficient_t;
				components<coefficient_t, PossibleGrades> result;
				if (lhs.basis_blade().get() == rhs.basis_blade().get()) {
					auto result_coefficient = lhs.coefficient() + rhs.coefficient();
					if (result_coefficient != 0) {
						result.insert(lhs.basis_blade(), result_coefficient);
					}
				}
				else {
					if (lhs.coefficient() != 0) {
						result.insert(lhs.basis_blade(), lhs.coefficient());
					}
					if (rhs.coefficient() != 0) {
						result.insert(rhs.basis_blade(), rhs.coefficient());
					}
				}
				return result;
				/*/
				components<decltype(LeftCoefficientType() + RightCoefficientType()), PossibleGrades> result;
				return result;
				/**/
			}

			template<class LeftCoefficientType, default_bitset_t PossibleGrades, class RightCoefficientType>
			constexpr decltype(auto) binary_plus_element(components<LeftCoefficientType, PossibleGrades> const &lhs, component<RightCoefficientType, dbasis_blade<PossibleGrades> > const &rhs) {
				//TODO LAZY O de cima � o original
				/**
				typedef typename std::common_type<LeftCoefficientType, RightCoefficientType, decltype(LeftCoefficientType() + RightCoefficientType())>::type coefficient_t;

				components<coefficient_t, PossibleGrades> result;

				auto lhs_itr = lhs.begin(), lhs_end = lhs.end();
				bool rhs_not_used = true;
				while (lhs_itr != lhs_end && rhs_not_used) {
					if (lhs_itr->first.get() < rhs.basis_blade().get()) {
						result.insert(lhs_itr->first, lhs_itr->second);
						++lhs_itr;
					}
					else if (lhs_itr->first.get() > rhs.basis_blade().get()) {
						if (rhs.coefficient() != 0) {
							result.insert(rhs.basis_blade(), rhs.coefficient());
						}
						rhs_not_used = false;
					}
					else {
						auto aux = lhs_itr->second + rhs.coefficient();
						if (aux != 0) {
							result.insert(lhs_itr->first, aux);
						}
						++lhs_itr;
						rhs_not_used = false;
					}
				}

				for (; lhs_itr != lhs_end; ++lhs_itr) {
					result.insert(lhs_itr->first, lhs_itr->second);
				}

				if (rhs_not_used && rhs.coefficient() != 0) {
					result.insert(rhs.basis_blade(), rhs.coefficient());
				}

				return result;
				/*/
				components<decltype(LeftCoefficientType() + RightCoefficientType()), PossibleGrades> result;
				return result;
				/**/
			}

			template<class LeftCoefficientType, default_bitset_t PossibleGrades, class RightCoefficientType>
			constexpr decltype(auto) binary_plus_element(component<LeftCoefficientType, dbasis_blade<PossibleGrades> > const &lhs, components<RightCoefficientType, PossibleGrades> const &rhs) {
				//TODO lazy
				typedef typename std::common_type<LeftCoefficientType, RightCoefficientType, decltype(LeftCoefficientType() + RightCoefficientType())>::type coefficient_t;

				components<coefficient_t, PossibleGrades> result;

				bool lhs_not_used = true;
				auto rhs_itr = rhs.begin(), rhs_end = rhs.end();
				while (lhs_not_used && rhs_itr != rhs_end) {
					if (lhs.basis_blade().get() < rhs_itr->first.get()) {
						if (lhs.coefficient() != 0) {
							result.insert(lhs.basis_blade(), lhs.coefficient());
						}
						lhs_not_used = false;
					}
					else if (lhs.basis_blade().get() > rhs_itr->first.get()) {
						result.insert(rhs_itr->first, rhs_itr->second);
						++rhs_itr;
					}
					else {
						auto aux = lhs.coefficient() + rhs_itr->second;
						if (aux != 0) {
							result.insert(lhs.basis_blade(), aux);
						}
						lhs_not_used = false;
						++rhs_itr;
					}
				}

				if (lhs_not_used && lhs.coefficient() != 0) {
					result.insert(lhs.basis_blade(), lhs.coefficient());
				}

				for (; rhs_itr != rhs_end; ++rhs_itr) {
					result.insert(rhs_itr->first, rhs_itr->second);
				}

				return result;
			}

			template<class LeftCoefficientType, default_bitset_t PossibleGrades, class RightCoefficientType>
			constexpr decltype(auto) binary_plus_element(components<LeftCoefficientType, PossibleGrades> const &lhs, components<RightCoefficientType, PossibleGrades> const &rhs) {
				//TODO lazy
				typedef typename std::common_type<LeftCoefficientType, RightCoefficientType, decltype(LeftCoefficientType() + RightCoefficientType())>::type coefficient_t;

				components<coefficient_t, PossibleGrades> result;

				auto lhs_itr = lhs.begin(), lhs_end = lhs.end();
				auto rhs_itr = rhs.begin(), rhs_end = rhs.end();
				while (lhs_itr != lhs_end && rhs_itr != rhs_end) {
					if (lhs_itr->first.get() < rhs_itr->first.get()) {
						result.insert(lhs_itr->first, lhs_itr->second);
						++lhs_itr;
					}
					else if (lhs_itr->first.get() > rhs_itr->first.get()) {
						result.insert(rhs_itr->first, rhs_itr->second);
						++rhs_itr;
					}
					else {
						auto aux = lhs_itr->second + rhs_itr->second;
						if (aux != 0) {
							result.insert(lhs_itr->first, aux);
						}
						++lhs_itr;
						++rhs_itr;
					}
				}

				for (; lhs_itr != lhs_end; ++lhs_itr) {
					result.insert(lhs_itr->first, lhs_itr->second);
				}

				for (; rhs_itr != rhs_end; ++rhs_itr) {
					result.insert(rhs_itr->first, rhs_itr->second);
				}

				return result;
			}

		}

	}

}

#endif // __GA_CLIFFORD_ELEMENT_BINARY_PLUS_HPP__