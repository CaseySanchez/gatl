_GA_DEFINE_BINARY_METRIC_OPERATION(gp)

template<class FirstArgumentType, class SecondArgumentType>
inline
auto operator * (FirstArgumentType&& m1, SecondArgumentType&& m2)->decltype(gp_em(m1, m2)) {
	///TODO Verificar se um dos dois argumentos � um escalar. Caso contr�rio, levantar ga::illegal_call_exception().
	return gp_em(m1, m2);
}
