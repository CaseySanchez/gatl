#ifndef __GA_CORE_LAZY_CONTEXT_HPP__
#define __GA_CORE_LAZY_CONTEXT_HPP__

namespace ga {

	namespace detail {
	
		// Returns the greater ID found in the given expressions.
		template<class... Expressions>
		struct greater_id;

		template<class... Expressions>
		constexpr tag_t greater_id_v = greater_id<Expressions...>::value;

		template<class Expression, class... NextExpressions>
		struct greater_id<Expression, NextExpressions...> {
			constexpr static tag_t value = greater(greater_id_v<Expression>, greater_id_v<NextExpressions...>); // recursion
		};

		template<>
		struct greater_id<> {
			constexpr static tag_t value = 0; // end of recursion
		};

		template<class Expression>
		struct greater_id<Expression> {
			constexpr static tag_t value = 0; // default
		};

		template<tag_t Tag, std::size_t Index>
		struct greater_id<get_value<Tag, Index> > {
			constexpr static tag_t value = Tag;
		};

		template<tag_t Tag, std::size_t Index>
		struct greater_id<get_map_values<Tag, Index> > {
			constexpr static tag_t value = Tag;
		};

		template<tag_t Tag, std::size_t Index>
		struct greater_id<get_bitset<Tag, Index> > {
			constexpr static tag_t value = Tag;
		};

		template<tag_t Tag, std::size_t Index>
		struct greater_id<get_map_bitsets<Tag, Index> > {
			constexpr static tag_t value = Tag;
		};

		template<default_bitset_t PossibleGrades, class Bitset>
		struct greater_id<dynamic_basis_blade<PossibleGrades, Bitset> > {
			constexpr static tag_t value = greater_id_v<Bitset>;
		};

		template<class Coefficient, class BasisBlade>
		struct greater_id<component<Coefficient, BasisBlade> > {
			constexpr static tag_t value = greater(greater_id_v<Coefficient>, greater_id_v<BasisBlade>);
		};

		template<name_t Name, class... Arguments>
		struct greater_id<function<Name, Arguments...> > {
			constexpr static tag_t value = greater_id_v<Arguments...>;
		};

		// Produces an expression where stored_value, stored_map_values, stored_bitset, and stored_map_bitsets change to get_value<Tag, Index>, get_map_values<Tag, Index>, get_bitset<Tag, Index>, and get_map_bitsets<Tag, Index>, respectively.
		template<class Expression, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables;

		template<class Expression, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		using tag_variables_t = typename tag_variables<Expression, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex>::type;

		template<class Expression, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables {
			typedef Expression type; // default
		};

		template<tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<stored_value, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef get_value<Tag, BaseValueIndex> type;
		};

		template<tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<stored_bitset, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef get_bitset<Tag, BaseBitsetIndex> type;
		};

		template<default_bitset_t PossibleGrades, class Bitset, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<dynamic_basis_blade<PossibleGrades, Bitset>, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef dynamic_basis_blade_t<PossibleGrades, tag_variables_t<Bitset, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> > type;
		};

		template<class Coefficient, class BasisBlade, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<component<Coefficient, BasisBlade>, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef component_t<tag_variables_t<Coefficient, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex>, tag_variables_t<BasisBlade, Tag, BaseValueIndex + count_stored_values_v<Coefficient>, BaseBitsetIndex + count_stored_bitsets_v<Coefficient>, BaseMapIndex + count_stored_maps_v<Coefficient> > > type;
		};

		template<default_bitset_t PossibleGrades, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<component<stored_map_values, dynamic_basis_blade<PossibleGrades, stored_map_bitsets> >, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef component_t<get_map_values<Tag, BaseMapIndex>, dynamic_basis_blade<PossibleGrades, get_map_bitsets<Tag, BaseMapIndex> > > type;
		};

		template<tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex, class Function, class... Arguments>
		struct _tag_variables_in_function;

		template<tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex, name_t Name, class... FunctionArguments, class Argument, class... NextArguments>
		struct _tag_variables_in_function<Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex, function<Name, FunctionArguments...>, Argument, NextArguments...> {
			typedef typename _tag_variables_in_function<Tag, BaseValueIndex + count_stored_values_v<Argument>, BaseBitsetIndex + count_stored_bitsets_v<Argument>, BaseMapIndex + count_stored_maps_v<Argument>, function<Name, FunctionArguments..., tag_variables_t<Argument, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> >, NextArguments...>::type type;
		};

		template<tag_t Tag, std::size_t ValueIndexEnd, std::size_t BitsetIndexEnd, std::size_t MapIndexEnd, name_t Name, class... FunctionArguments>
		struct _tag_variables_in_function<Tag, ValueIndexEnd, BitsetIndexEnd, MapIndexEnd, function<Name, FunctionArguments...> > {
			typedef function<Name, FunctionArguments...> type;
		};

		template<name_t Name, class... Arguments, tag_t Tag, std::size_t BaseValueIndex, std::size_t BaseBitsetIndex, std::size_t BaseMapIndex>
		struct tag_variables<function<Name, Arguments...>, Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex> {
			typedef typename _tag_variables_in_function<Tag, BaseValueIndex, BaseBitsetIndex, BaseMapIndex, function<Name>, Arguments...>::type type;
		};

		// Evaluates the given (or the collection of) clifford_expression<...>.
		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct eval_clifford_expression;

		template<tag_t LowerTag, tag_t UpperTag, class Expression, class... InputTypes>
		using eval_coefficient_t = typename eval_clifford_expression<LowerTag, UpperTag, Expression>::template coefficient_type<InputTypes...>::type;

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		using eval_expression_t = typename eval_clifford_expression<LowerTag, UpperTag, Expression>::expression_type;

		template<tag_t LowerTag, tag_t UpperTag, class Expression, class... NextExpressions>
		struct eval_clifford_expressions {
			template<class... InputTypes>
			struct coefficient_type {
				typedef std::common_type_t<typename eval_clifford_expression<LowerTag, UpperTag, Expression>::template coefficient_type<InputTypes...>::type, typename eval_clifford_expressions<LowerTag, UpperTag, NextExpressions...>::template coefficient_type<InputTypes...>::type> type;
			};

			// expression_type is not defined here.

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr &value_itr, BitsetItr &bitset_itr, MapIts &map_itr, std::tuple<InputTypes...> const &args) {
				eval_clifford_expression<LowerTag, UpperTag, Expression>::run(value_itr, bitset_itr, map_itr, args);
				eval_clifford_expressions<LowerTag, UpperTag, NextExpressions...>::run(value_itr, bitset_itr, map_itr, args);
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct eval_clifford_expressions<LowerTag, UpperTag, Expression> {
			template<class... InputTypes>
			struct coefficient_type {
				typedef typename eval_clifford_expression<LowerTag, UpperTag, Expression>::template coefficient_type<InputTypes...>::type type;
			};

			// expression_type is not defined here.

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr &value_itr, BitsetItr &bitset_itr, MapIts &map_itr, std::tuple<InputTypes...> const &args) {
				eval_clifford_expression<LowerTag, UpperTag, Expression>::run(value_itr, bitset_itr, map_itr, args);
			}
		};

		template<class Expression>
		struct _eval_clifford_expression_do_nothing {
			template<class... InputTypes>
			struct coefficient_type {
				typedef default_integral_t type;
			};

			typedef Expression expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr const &, BitsetItr const &, MapIts const &, std::tuple<InputTypes...> const &) {
				// Do nothing.
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class ExpressionType, class... Expressions>
		struct _eval_clifford_expression_move {
			template<class... InputTypes>
			struct coefficient_type {
				typedef typename eval_clifford_expressions<LowerTag, UpperTag, Expressions...>::template coefficient_type<InputTypes...>::type type;
			};

			typedef ExpressionType expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr &value_itr, BitsetItr &bitset_itr, MapIts &map_itr, std::tuple<InputTypes...> const &args) {
				eval_clifford_expressions<LowerTag, UpperTag, Expressions...>::run(value_itr, bitset_itr, map_itr, args);
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct _eval_clifford_expression_store_value {
			template<class... InputTypes>
			struct coefficient_type {
				typedef decltype(Expression::eval<LowerTag, UpperTag>(std::declval<std::tuple<InputTypes...> >())) type;
			};

			typedef stored_value expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr &value_itr, BitsetItr const &, MapIts const &, std::tuple<InputTypes...> const &args) {
				*value_itr = Expression::eval<LowerTag, UpperTag>(args);
				std::advance(value_itr, 1);
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct _eval_clifford_expression_store_map_values {
			template<class... InputTypes>
			struct coefficient_type {
				typedef nullptr_t type; //TODO Not supported yet (map)
			};

			typedef stored_map_values expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr const &, BitsetItr const &, MapIts &map_itr, std::tuple<InputTypes...> const &args) {
				//TODO Not supported yet (map)
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct _eval_clifford_expression_store_bitset {
			template<class... InputTypes>
			struct coefficient_type {
				typedef default_integral_t type;
			};

			typedef stored_bitset expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr const &, BitsetItr &bitset_itr, MapIts const &, std::tuple<InputTypes...> const &args) {
				*bitset_itr = Expression::eval<LowerTag, UpperTag>(args);
				std::advance(bitset_itr, 1);
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct _eval_clifford_expression_store_map_bitsets {
			template<class... InputTypes>
			struct coefficient_type {
				typedef nullptr_t type; //TODO Not supported yet (map)
			};

			typedef stored_map_bitsets expression_type;

			template<class ValueItr, class BitsetItr, class MapIts, class... InputTypes>
			constexpr static void run(ValueItr const &, BitsetItr const &, MapIts &map_itr, std::tuple<InputTypes...> const &args) {
				//TODO Not supported yet (map)
			}
		};

		template<tag_t LowerTag, tag_t UpperTag, class Expression>
		struct eval_clifford_expression : _eval_clifford_expression_do_nothing<Expression> {
		};

		template<tag_t LowerTag, tag_t UpperTag, tag_t Tag, std::size_t Index>
		struct eval_clifford_expression<LowerTag, UpperTag, get_value<Tag, Index> > :
			std::conditional_t<
				(LowerTag <= Tag && Tag <= UpperTag),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, get_value<Tag, Index> >,
				_eval_clifford_expression_do_nothing<get_value<Tag, Index> >
			> {
		};

		template<tag_t LowerTag, tag_t UpperTag, tag_t Tag, std::size_t Index>
		struct eval_clifford_expression<LowerTag, UpperTag, get_map_values<Tag, Index> > :
			std::conditional_t<
				(LowerTag <= Tag && Tag <= UpperTag),
				_eval_clifford_expression_store_map_values<LowerTag, UpperTag, get_map_values<Tag, Index> >,
				_eval_clifford_expression_do_nothing<get_map_values<Tag, Index> >
			> {
		};

		template<tag_t LowerTag, tag_t UpperTag, tag_t Tag, std::size_t Index>
		struct eval_clifford_expression<LowerTag, UpperTag, get_bitset<Tag, Index> > :
			std::conditional_t<
				(LowerTag <= Tag && Tag <= UpperTag),
				_eval_clifford_expression_store_bitset<LowerTag, UpperTag, get_bitset<Tag, Index> >,
				_eval_clifford_expression_do_nothing<get_bitset<Tag, Index> >
			> {
		};

		template<tag_t LowerTag, tag_t UpperTag, tag_t Tag, std::size_t Index>
		struct eval_clifford_expression<LowerTag, UpperTag, get_map_bitsets<Tag, Index> > :
			std::conditional_t<
				(LowerTag <= Tag && Tag <= UpperTag),
				_eval_clifford_expression_store_map_bitsets<LowerTag, UpperTag, get_map_bitsets<Tag, Index> >,
				_eval_clifford_expression_do_nothing<get_map_bitsets<Tag, Index> >
			> {
		};

		template<tag_t LowerTag, tag_t UpperTag, default_bitset_t PossibleGrades, class Bitset>
		struct eval_clifford_expression<LowerTag, UpperTag, dynamic_basis_blade<PossibleGrades, Bitset> > :
			_eval_clifford_expression_move<LowerTag, UpperTag, dynamic_basis_blade_t<PossibleGrades, eval_expression_t<LowerTag, UpperTag, Bitset> >, Bitset> {
		};

		template<tag_t LowerTag, tag_t UpperTag, class Coefficient, class BasisBlade>
		struct eval_clifford_expression<LowerTag, UpperTag, component<Coefficient, BasisBlade> > :
			_eval_clifford_expression_move<LowerTag, UpperTag, component_t<eval_expression_t<LowerTag, UpperTag, Coefficient>, eval_expression_t<LowerTag, UpperTag, BasisBlade> >, Coefficient, BasisBlade> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, name_t Name, class Argument>
		struct eval_clifford_expression<LowerTag, UpperTag, function<Name, Argument> > :
			std::conditional_t<
				std::is_same_v<eval_expression_t<LowerTag, UpperTag, Argument>, stored_value>,
				_eval_clifford_expression_store_value<LowerTag, UpperTag, function<Name, Argument> >,
				std::conditional_t<
					std::is_same_v<eval_expression_t<LowerTag, UpperTag, Argument>, stored_bitset>,
					_eval_clifford_expression_store_bitset<LowerTag, UpperTag, function<Name, Argument> >,
					_eval_clifford_expression_move<LowerTag, UpperTag, typename function<Name, eval_expression_t<LowerTag, UpperTag, Argument> >::type, Argument>
				>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, name_t Name, class LeftArgument, class RightArgument>
		struct eval_clifford_expression<LowerTag, UpperTag, function<Name, LeftArgument, RightArgument> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, LeftArgument>, stored_value> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, RightArgument> >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, LeftArgument> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, RightArgument>, stored_value>),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, function<Name, LeftArgument, RightArgument> >,
				std::conditional_t<
					(std::is_same_v<eval_expression_t<LowerTag, UpperTag, LeftArgument>, stored_bitset> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, RightArgument> >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, LeftArgument> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, RightArgument>, stored_bitset>),
					_eval_clifford_expression_store_bitset<LowerTag, UpperTag, function<Name, LeftArgument, RightArgument> >,
					_eval_clifford_expression_move<LowerTag, UpperTag, typename function<Name, eval_expression_t<LowerTag, UpperTag, LeftArgument>, eval_expression_t<LowerTag, UpperTag, RightArgument> >::type, LeftArgument, RightArgument>
				>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Argument, class... NextArguments>
		struct eval_clifford_expression<LowerTag, UpperTag, add<Argument, NextArguments...> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, Argument>, stored_value> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, add_t<NextArguments...> > >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, Argument> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, add_t<NextArguments...> >, stored_value>),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, add<Argument, NextArguments...> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, addition_t<eval_expression_t<LowerTag, UpperTag, Argument>, eval_expression_t<LowerTag, UpperTag, add_t<NextArguments...> > >, Argument, add_t<NextArguments...> >
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Coefficient, class BasisBlade, class... NextComponents>
		struct eval_clifford_expression<LowerTag, UpperTag, add<component<Coefficient, BasisBlade>, NextComponents...> > :
			_eval_clifford_expression_move<LowerTag, UpperTag, addition_t<component_t<eval_expression_t<LowerTag, UpperTag, Coefficient>, eval_expression_t<LowerTag, UpperTag, BasisBlade> >, eval_expression_t<LowerTag, UpperTag, add_t<NextComponents...> > >, Coefficient, BasisBlade, add_t<NextComponents...> > {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Argument, class... NextArguments>
		struct eval_clifford_expression<LowerTag, UpperTag, mul<Argument, NextArguments...> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, Argument>, stored_value> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, mul_t<NextArguments...> > >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, Argument> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, mul_t<NextArguments...> >, stored_value>),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, mul<Argument, NextArguments...> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, product_t<eval_expression_t<LowerTag, UpperTag, Argument>, eval_expression_t<LowerTag, UpperTag, mul_t<NextArguments...> >, value_mapping>, Argument, mul_t<NextArguments...> >
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class LeftArgument, class RightArgument>
		struct eval_clifford_expression<LowerTag, UpperTag, power<LeftArgument, RightArgument> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, LeftArgument>, stored_value> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, RightArgument> >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, LeftArgument> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, RightArgument>, stored_value>),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, power<LeftArgument, RightArgument> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, power_t<eval_expression_t<LowerTag, UpperTag, LeftArgument>, eval_expression_t<LowerTag, UpperTag, RightArgument> >, LeftArgument, RightArgument>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class LeftBitset, class RightBitset>
		struct eval_clifford_expression<LowerTag, UpperTag, reordering_sign<LeftBitset, RightBitset> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, LeftBitset>, stored_bitset> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, RightBitset> >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, LeftBitset> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, RightBitset>, stored_bitset>),
				_eval_clifford_expression_store_value<LowerTag, UpperTag, reordering_sign<LeftBitset, RightBitset> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, reordering_sign_t<eval_expression_t<LowerTag, UpperTag, LeftBitset>, eval_expression_t<LowerTag, UpperTag, RightBitset> >, LeftBitset, RightBitset>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Bitset>
		struct eval_clifford_expression<LowerTag, UpperTag, count_one_bits<Bitset> > :
			std::conditional_t<
				std::is_same_v<eval_expression_t<LowerTag, UpperTag, Bitset>, stored_bitset>,
				_eval_clifford_expression_store_value<LowerTag, UpperTag, count_one_bits<Bitset> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, count_one_bits_t<eval_expression_t<LowerTag, UpperTag, Bitset> >, Bitset>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class LeftBitset, class RightValue>
		struct eval_clifford_expression<LowerTag, UpperTag, bitwise_left_shift<LeftBitset, RightValue> > :
			std::conditional_t<
				(std::is_same_v<eval_expression_t<LowerTag, UpperTag, LeftBitset>, stored_bitset> && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, RightValue> >) || (can_be_stored_v<eval_expression_t<LowerTag, UpperTag, LeftBitset> > && std::is_same_v<eval_expression_t<LowerTag, UpperTag, RightValue>, stored_value>),
				_eval_clifford_expression_store_bitset<LowerTag, UpperTag, bitwise_left_shift<LeftBitset, RightValue> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, bitwise_left_shift_t<eval_expression_t<LowerTag, UpperTag, LeftBitset>, eval_expression_t<LowerTag, UpperTag, RightValue> >, LeftBitset, RightValue>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Bitset>
		struct eval_clifford_expression<LowerTag, UpperTag, bitwise_uminus<Bitset> > :
			std::conditional_t<
				std::is_same_v<eval_expression_t<LowerTag, UpperTag, Bitset>, stored_bitset>,
				_eval_clifford_expression_store_value<LowerTag, UpperTag, bitwise_uminus<Bitset> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, bitwise_uminus_t<eval_expression_t<LowerTag, UpperTag, Bitset> >, Bitset>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Bitset>
		struct eval_clifford_expression<LowerTag, UpperTag, bitwise_dec<Bitset> > :
			std::conditional_t<
				std::is_same_v<eval_expression_t<LowerTag, UpperTag, Bitset>, stored_bitset>,
				_eval_clifford_expression_store_value<LowerTag, UpperTag, bitwise_dec<Bitset> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, bitwise_dec_t<eval_expression_t<LowerTag, UpperTag, Bitset> >, Bitset>
			> {
		};

		//TODO Not supported yet (map)
		template<tag_t LowerTag, tag_t UpperTag, class Test, class TrueValue, class FalseValue>
		struct eval_clifford_expression<LowerTag, UpperTag, if_else<Test, TrueValue, FalseValue> > :
			std::conditional_t<
				can_be_stored_v<eval_expression_t<LowerTag, UpperTag, Test> > && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, TrueValue> > && can_be_stored_v<eval_expression_t<LowerTag, UpperTag, TrueValue> >,
				_eval_clifford_expression_store_value<LowerTag, UpperTag, if_else<Test, TrueValue, FalseValue> >,
				_eval_clifford_expression_move<LowerTag, UpperTag, if_else_t<eval_expression_t<LowerTag, UpperTag, Test>, eval_expression_t<LowerTag, UpperTag, TrueValue>, eval_expression_t<LowerTag, UpperTag, TrueValue> >, Test, TrueValue, FalseValue>
			> {
		};
		
		template<tag_t LowerTag, tag_t UpperTag, class CoefficientType, class Expression, class... InputTypes>
		constexpr static decltype(auto) eval(clifford_expression<CoefficientType, Expression> const &expression, std::tuple<InputTypes...> const &args) {
			typedef clifford_expression<eval_coefficient_t<LowerTag, UpperTag, Expression, std::remove_const_t<std::remove_reference_t<InputTypes> >...>, eval_expression_t<LowerTag, UpperTag, Expression> > result_type;

			typename result_type::value_storage_type values;
			typename result_type::bitset_storage_type bitsets;
			typename result_type::map_storage_type maps;

			auto value_itr = values.begin();
			auto bitset_itr = bitsets.begin();
			auto map_itr = maps.begin();

			eval_clifford_expression<LowerTag, UpperTag, Expression>::run(value_itr, bitset_itr, map_itr, args);

			return result_type(std::move(values), std::move(bitsets), std::move(maps));
		}

		// Superclass for ga::lazy_context<InputTypes...>.
		template<tag_t BaseTag, class... InputTypes>
		class _super_lazy_context;

		template<size_t ReverseIndex, class InputCoefficientType, class InputExpression, bool StoredReference = has_stored_entries_v<InputExpression> >
		class _super_lazy_context_input {
		public:

			typedef clifford_expression<InputCoefficientType, InputExpression> input_type;

			constexpr _super_lazy_context_input(_super_lazy_context_input const &) = default;
			constexpr _super_lazy_context_input(_super_lazy_context_input &&) = default;

			constexpr _super_lazy_context_input(input_type const &input) :
				input_(input) {
			}

			constexpr _super_lazy_context_input & operator=(_super_lazy_context_input const &) = delete;
			constexpr _super_lazy_context_input & operator=(_super_lazy_context_input &&) = delete;

			constexpr decltype(auto) get_as_tuple() const {
				return std::tie(input_);
			}

			constexpr static bool is_stored() {
				return true;
			}

		private:

			input_type const &input_;
		};

		template<std::size_t ReverseIndex, class InputCoefficientType, class InputExpression>
		class _super_lazy_context_input<ReverseIndex, InputCoefficientType, InputExpression, false> {
		public:

			typedef clifford_expression<InputCoefficientType, InputExpression> input_type;

			constexpr _super_lazy_context_input(_super_lazy_context_input const &) = default;
			constexpr _super_lazy_context_input(_super_lazy_context_input &&) = default;

			constexpr _super_lazy_context_input(clifford_expression<InputCoefficientType, InputExpression> const &) {
			}

			constexpr _super_lazy_context_input & operator=(_super_lazy_context_input const &) = delete;
			constexpr _super_lazy_context_input & operator=(_super_lazy_context_input &&) = delete;

			constexpr decltype(auto) get_as_tuple() const {
				return std::make_tuple();
			}

			constexpr static bool is_stored() {
				return false;
			}
		};

		template<tag_t BaseTag, class InputCoefficientType, class InputExpression, class... OtherInputCoefficientTypes, class... OtherInputExpressions>
		class _super_lazy_context<BaseTag, clifford_expression<InputCoefficientType, InputExpression>, clifford_expression<OtherInputCoefficientTypes, OtherInputExpressions>...> :
			private _super_lazy_context_input<sizeof...(OtherInputExpressions), InputCoefficientType, InputExpression>,
			private _super_lazy_context<_super_lazy_context_input<sizeof...(OtherInputExpressions), InputCoefficientType, InputExpression>::is_stored() ? BaseTag + 1 : BaseTag, clifford_expression<OtherInputCoefficientTypes, OtherInputExpressions>...> {
		private:

			typedef _super_lazy_context_input<sizeof...(OtherInputExpressions), InputCoefficientType, InputExpression> super_input;
			typedef _super_lazy_context<super_input::is_stored() ? BaseTag + 1 : BaseTag, clifford_expression<OtherInputCoefficientTypes, OtherInputExpressions>...> super_recursive;

		public:

			template<std::size_t Index>
			struct argument {
				typedef typename super_recursive::template argument<Index - 1>::type type;
			};

			template<>
			struct argument<0> {
				typedef std::conditional_t<
					super_input::is_stored(),
					clifford_expression<InputCoefficientType, tag_variables_t<InputExpression, BaseTag + 1, 0, 0, 0> >,
					typename super_input::input_type
				> type;
			};

			constexpr _super_lazy_context(_super_lazy_context const &) = default;
			constexpr _super_lazy_context(_super_lazy_context &&) = default;

			constexpr _super_lazy_context(clifford_expression<InputCoefficientType, InputExpression> const &input, clifford_expression<OtherInputCoefficientTypes, OtherInputExpressions> const &... other_inputs) :
				super_input(input),
				super_recursive(other_inputs...) {
			}

			constexpr _super_lazy_context & operator=(_super_lazy_context const &) = delete;
			constexpr _super_lazy_context & operator=(_super_lazy_context &&) = delete;

			constexpr decltype(auto) stored_inputs_tuple() const {
				return std::tuple_cat(super_input::get_as_tuple(), super_recursive::stored_inputs_tuple());
			}

			constexpr static std::size_t stored_inputs_count() {
				return (super_input::is_stored() ? 1 : 0) + super_recursive::stored_inputs_count();
			}
		};

		template<tag_t BaseTag>
		class _super_lazy_context<BaseTag> {
		public:

			constexpr _super_lazy_context() = default;
			constexpr _super_lazy_context(_super_lazy_context const &) = default;
			constexpr _super_lazy_context(_super_lazy_context &&) = default;

			constexpr _super_lazy_context & operator=(_super_lazy_context const &) = delete;
			constexpr _super_lazy_context & operator=(_super_lazy_context &&) = delete;

			constexpr decltype(auto) stored_inputs_tuple() const {
				return std::make_tuple();
			}

			constexpr static std::size_t stored_inputs_count() {
				return 0;
			}
		};
	
	}

	// Helper structure to define lazy arguments for lazy evaluation of expressions.
	template<class... InputTypes>
	class lazy_context;

	template<class... InputCoefficientTypes, class... InputExpressions>
	class lazy_context<clifford_expression<InputCoefficientTypes, InputExpressions>...> final :
		private detail::_super_lazy_context<detail::greater_id_v<InputExpressions...>, clifford_expression<InputCoefficientTypes, InputExpressions>...> {
	private:

		constexpr static detail::tag_t base_id = detail::greater_id_v<InputExpressions...>;
		
		typedef detail::_super_lazy_context<base_id, clifford_expression<InputCoefficientTypes, InputExpressions>...> super;

	public:

		template<std::size_t Index>
		using argument_t = typename super::template argument<Index>::type;

		template<std::size_t Index>
		using argument_coefficient_t = typename super::template argument<Index>::type::coefficient_type;

		template<std::size_t Index>
		using argument_expression_t = typename super::template argument<Index>::type::expression_type;

		constexpr lazy_context(lazy_context const &) = default;
		constexpr lazy_context(lazy_context &&) = default;

		constexpr lazy_context(clifford_expression<InputCoefficientTypes, InputExpressions> const &... inputs) :
			super(inputs...) {
		}

		constexpr lazy_context & operator=(lazy_context const &) = delete;
		constexpr lazy_context & operator=(lazy_context &&) = delete;

		template<std::size_t Index>
		constexpr static decltype(auto) argument() {
			return argument_t<Index>();
		}

		template<class CoefficientType, class Expression, class = std::enable_if_t<super::stored_inputs_count() != 0> >
		constexpr decltype(auto) eval(clifford_expression<CoefficientType, Expression> const &expression) const {
			return detail::eval<base_id + 1, base_id + (detail::tag_t)super::stored_inputs_count()>(expression, super::stored_inputs_tuple());
		}

		template<class CoefficientType, class Expression, class = std::enable_if_t<super::stored_inputs_count() == 0> >
		constexpr decltype(auto) eval(clifford_expression<CoefficientType, Expression> &&) const {
			return clifford_expression<CoefficientType, Expression>();
		}
	};

	template<class... InputTypes>
	constexpr decltype(auto) make_lazy_context(InputTypes const &... inputs) {
		return lazy_context<InputTypes...>(inputs...);
	}

}

#endif // __GA_CORE_LAZY_CONTEXT_HPP__
