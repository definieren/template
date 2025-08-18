```cpp
struct Combination {
	int N;
	std::vector<Z> _fac, _ifac, _inv;
	
	Combination(int n): N(0), _fac{1}, _ifac{1}, _inv{0} { init(n); }
	Combination(): N(0), _fac{1}, _ifac{1}, _inv{0} {}
	void init(int n) {
		if (n <= N) return;
		_fac.resize(n + 1), _ifac.resize(n + 1), _inv.resize(n + 1);
		for (int i = N + 1; i <= n; i ++) _fac[i] = _fac[i - 1] * i;
		_ifac[n] = _fac[n].inv();
		for (int i = n; i > N; i --) _ifac[i - 1] = _ifac[i] * i,
																 _inv[i] = _ifac[i] * _fac[i - 1];
		N = n; return;
	}
	Z fac(int n) {
		init(n << 1);
		return n < 0 ? Z::from_raw(0) : _fac[n];
	}
	Z ifac(int n) {
		init(n << 1);
		return n < 0 ? Z::from_raw(0) : _ifac[n];
	}
	Z inv(int n) {
		init(n << 1);
		return n < 0 ? Z::from_raw(0) : _inv[n];
	}
	Z binom(int n, int m) {
		if (n < m || n < 0 || m < 0) return 0;
		return fac(n) * ifac(m) * ifac(n - m);
	}
} comb;
```