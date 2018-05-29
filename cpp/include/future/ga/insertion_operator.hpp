#ifndef __FUTURE_GA_INSERTION_OPERATOR_HPP__
#define __FUTURE_GA_INSERTION_OPERATOR_HPP__

namespace ga {

	namespace detail {

		void write_basis_blade(std::ostream &os, default_bitset_t arg) {
			if (arg == default_bitset_t(0)) {
				os << "1";
			}
			else {
				index_t ind = 0;
				bool first = true;
				while (arg != default_bitset_t(0)) {
					if ((arg & default_bitset_t(1)) != default_bitset_t(0)) {
						if (!first) os << "^";
						else first = false;
						os << "e" << (ind + 1);
					}
					arg >>= 1;
					++ind;
				}
			}
		}

		template<class Expression>
		struct write_expression;

		template<default_integral_t Value>
		struct write_expression<constant_value<Value> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "<" << Value << ">";
				first = false;
			}
		};

		template<tag_t Tag, std::size_t Index>
		struct write_expression<get_value<Tag, Index> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "{Tag: " << Tag << ", ValueIndex: " << Index << "}";
				first = false;
			}
		};

		template<tag_t Tag, std::size_t Index>
		struct write_expression<get_map_values<Tag, Index> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "[{Tag: " << Tag << ", MapValuesIndex: " << Index << "}]";
				first = false;
			}
		};

		template<>
		struct write_expression<stored_value> {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &, MapCIts &, bool &first) {
				if (first) {
					os << (*value_citr);
					first = false;
				}
				else {
					if ((*value_citr) >= 0) os << (*value_citr);
					else os << "(" << (*value_citr) << ")";
				}
				std::advance(value_citr, 1);
			}
		};

		template<default_bitset_t Bitset>
		struct write_expression<constant_bitset<Bitset> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "<";
				write_basis_blade(os, Bitset);
				os << ">";
				first = false;
			}
		};

		template<tag_t Tag, std::size_t Index>
		struct write_expression<get_bitset<Tag, Index> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "{Tag: " << Tag << ", BitsetIndex: " << Index << "}";
				first = false;
			}
		};

		template<tag_t Tag, std::size_t Index>
		struct write_expression<get_map_bitsets<Tag, Index> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &, bool &first) {
				os << "[{Tag: " << Tag << ", MapBitsetsIndex: " << Index << "}]";
				first = false;
			}
		};

		template<>
		struct write_expression<stored_bitset> {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &bitset_citr, MapCIts &, bool &first) {
				write_basis_blade(os, *bitset_citr);
				std::advance(bitset_citr, 1);
				first = false;
			}
		};

		template<default_bitset_t BasisVectors>
		struct write_expression<constant_basis_blade<BasisVectors> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				write_expression<constant_bitset<BasisVectors> >::run(os, value_citr, bitset_citr, map_citr, first);
				first = false;
			}
		};

		template<default_bitset_t PossibleGrades, class Bitset>
		struct write_expression<dynamic_basis_blade<PossibleGrades, Bitset> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				write_expression<Bitset>::run(os, value_citr, bitset_citr, map_citr, first);
				first = false;
			}
		};

		template<class Coefficient, class BasisBlade>
		struct write_expression<component<Coefficient, BasisBlade> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<Coefficient>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ") * ";
				write_expression<BasisBlade>::run(os, value_citr, bitset_citr, map_citr, first);
				first = false;
			}
		};

		template<default_bitset_t PossibleGrades>
		struct write_expression<component<stored_map_values, dynamic_basis_blade<PossibleGrades, stored_map_bitsets> > > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &, BitsetCItr &, MapCIts &map_citr, bool &first) {
				os << "[";
				if (!map_citr->empty()) {
					bool local_first = true;
					for (auto &curr : *map_citr) {
						if (!local_first) os << " + ";
						else local_first = false;

						if (curr.second >= 0) os << curr.second;
						else os << "(" << curr.second << ")";

						os << " * ";
						write_basis_blade(os, curr.first);
					}
				}
				else {
					os << c<0>;
				}
				os << "]";
				std::advance(map_citr, 1);
				first = false;
			}
		};

		template<class Argument, class... NextArguments>
		struct write_expression<add<Argument, NextArguments...> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				if (!first) os << " + ";
				write_expression<Argument>::run(os, value_citr, bitset_citr, map_citr, first);
				write_expression<add<NextArguments...> >::run(os, value_citr, bitset_citr, map_citr, first);
			}
		};

		template<class LeftArgument, class RightArgument>
		struct write_expression<add<LeftArgument, RightArgument> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				if (!first) os << " + ";
				write_expression<LeftArgument>::run(os, value_citr, bitset_citr, map_citr, first);
				os << " + ";
				write_expression<RightArgument>::run(os, value_citr, bitset_citr, map_citr, first);
			}
		};

		template<class Argument, class... NextArguments>
		struct write_expression<mul<Argument, NextArguments...> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				if (!first) os << " * ";
				write_expression<Argument>::run(os, value_citr, bitset_citr, map_citr, first);
				write_expression<mul<NextArguments...> >::run(os, value_citr, bitset_citr, map_citr, first);
			}
		};

		template<class LeftArgument, class RightArgument>
		struct write_expression<mul<LeftArgument, RightArgument> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				if (!first) os << " * ";
				write_expression<LeftArgument>::run(os, value_citr, bitset_citr, map_citr, first);
				os << " * ";
				write_expression<RightArgument>::run(os, value_citr, bitset_citr, map_citr, first);
			}
		};

		template<class LeftArgument, class RightArgument>
		struct write_expression<power<LeftArgument, RightArgument> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "pow(";
				write_expression<LeftArgument>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ", ";
				local_first = true;
				write_expression<RightArgument>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class LeftBitset, class RightBitset>
		struct write_expression<reordering_sign<LeftBitset, RightBitset> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "reordering_sign(";
				write_expression<LeftBitset>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ", ";
				local_first = true;
				write_expression<RightBitset>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};
	
		template<class Bitset>
		struct write_expression<count_one_bits<Bitset> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "count_one_bits(";
				write_expression<Bitset>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class LeftType, class RightValue>
		struct write_expression<bitwise_left_shift<LeftType, RightValue> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<LeftType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " LSHIFTb ";
				local_first = true;
				write_expression<RightValue>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class LeftType, class RightType>
		struct write_expression<bitwise_and<LeftType, RightType> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<LeftType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " ANDb ";
				local_first = true;
				write_expression<RightType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class LeftType, class RightType>
		struct write_expression<bitwise_xor<LeftType, RightType> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<LeftType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " XORb ";
				local_first = true;
				write_expression<RightType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class LeftType, class RightType>
		struct write_expression<equal<LeftType, RightType> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<LeftType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " == ";
				local_first = true;
				write_expression<RightType>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

		template<class Test, class TrueValue, class FalseValue>
		struct write_expression<if_else<Test, TrueValue, FalseValue> > {
			template<class ValueCItr, class BitsetCItr, class MapCIts>
			inline static void run(std::ostream &os, ValueCItr &value_citr, BitsetCItr &bitset_citr, MapCIts &map_citr, bool &first) {
				bool local_first = true;
				os << "(";
				write_expression<Test>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " ? ";
				local_first = true;
				write_expression<TrueValue>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << " : ";
				local_first = true;
				write_expression<FalseValue>::run(os, value_citr, bitset_citr, map_citr, local_first);
				os << ")";
				first = false;
			}
		};

	}

	template<class RightCoefficientType, class RightExpression>
	std::ostream & operator<<(std::ostream &os, clifford_expression<RightCoefficientType, RightExpression> const &rhs) {
		bool first = true;
		
		auto value_citr = rhs.values().cbegin();
		auto bitset_citr = rhs.bitsets().cbegin();
		auto map_citr = rhs.maps().cbegin();
		detail::write_expression<RightExpression>::run(os, value_citr, bitset_citr, map_citr, first);

		if (first) {
			os << c<0>;
		}

		return os;
	}

}

#endif // __FUTURE_GA_INSERTION_OPERATOR_HPP__