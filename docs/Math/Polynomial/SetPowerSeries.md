```cpp
template<const unsigned short N = 21>
struct SetPowerSeries : public vector<Z> {
	private:
		static constexpr void Ranked_Zeta(const int &n, const Z *f, Z (*F)[N + 1]) {
			for (int S = 0; S < 1 << n; S ++)
				fill(F[S], F[S] + n + 1, 0),
				F[S][__builtin_popcount(S)] = f[S];
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						for (int k = 0; k <= n; k ++)
							F[S | i][k] += F[S][k];
			return;
		}
		static constexpr void Ranked_Mobius(const int &n, Z (*F)[N + 1], Z *f) {
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						for (int k = 0; k <= n; k ++)
							F[S | i][k] -= F[S][k];
			for (int S = 0; S < 1 << n; S ++)
				f[S] = F[S][__builtin_popcount(S)];
			return;
		}
		static constexpr void Convolution_Naive(const int &n, const Z *f, const Z *g, Z *ret) {
			for (int S = 0, T; S < 1 << n; S ++)
				for (ret[T = S] = f[S] * g[0]; T; (-- T) &= S)
					ret[S] += f[S ^ T] * g[T];
			return;
		}
		static void Convolution_FWT(const int &n, const Z *f, const Z *g, Z *ret) {
			static Z F[1 << N | 32][N + 1], G[1 << N | 32][N + 1];
			Ranked_Zeta(n, f, F), Ranked_Zeta(n, g, G);
			for (int S = 0; S < 1 << n; S ++) {
				Z x;
				for (int c = __builtin_popcount(S), i = min(2 * c, n), j; i >= c; F[S][i --] = x)
					for (x.x = 0, j = i - c; j <= c; j ++)
						x += F[S][j] * G[S][i - j];
			}
			Ranked_Mobius(n, F, ret);
			return;
		}
		static constexpr void Convolution(const int &n, const Z *f, const Z *g, Z *ret) {
			(n <= 10) ? Convolution_Naive(n, f, g, ret)
								: Convolution_FWT(n, f, g, ret);
			return;
		}
		static constexpr void Division_Naive(const int &n, Z *f, const Z *g) {
			for (int S = 1; S < 1 << n; S ++)
				for (int T = S; T; (-- T) &= S)
					f[S] -= f[S ^ T] * g[T];
			return;
		}
		static void Division_FWT(const int &n, Z *f, const Z *g) {
			static Z F[1 << N | 32][N + 1], G[1 << N | 32][N + 1];
			Ranked_Zeta(n, f, F), Ranked_Zeta(n, g, G);
			for (int S = 0; S < 1 << n; S ++)
				for (int c = __builtin_popcount(S), i = 0; i <= n; i ++)
					for (int j = max(0, i - c); j < i; j ++)
						F[S][i] -= F[S][j] * G[S][i - j];
			Ranked_Mobius(n, F, f);
			return;
		}
		static constexpr void Division(const int &n, Z *f, const Z *g) {
			(n <= 10) ? Division_Naive(n, f, g)
								: Division_FWT(n, f, g);
			return;
		}
		static constexpr void Composite(int n, const Z *egf, const Z *f, Z *g) {
			for (int i = n; ~i; i --) {
				for (int j = n - i; -- j >= 0; )
					Convolution(j, g, f + (1 << j), g + (1 << j));
				g[0] = egf[i];
			}
			return;
		}
		template<class Func = void(*)(int, Z&)>
		static constexpr void Online_Convolution_Naive(const int &n, Z *f, const Z *g, const Z &ori, const Func &func = [](int, Z&) {}) {
			f[0] = ori;
			for (int S = 1; S < 1 << n; func(S, f[S]), S ++)
				for (int T = S; T; (-- T) &= S)
					f[S] += f[S ^ T] * g[T];
			return;
		}
		template<class Func = void(*)(int, Z&)>
		static void Online_Convolution_FWT(const int &n, Z *f, const Z *g, const Z &ori, const Func &func = [](int, Z&) {}) {
			static Z F[N + 1][1 << N | 32], G[1 << N | 32][N + 1];
			fill(F[0], F[0] + (1 << n), ori), Ranked_Zeta(n, g, G);
			fill(F[1], F[1] + (1 << n), Z(0));
			for (int i = 1; i < 1 << n; i <<= 1)
				func(i, F[1][i] = ori * g[i]);
			for (int d = 2; d <= n; d ++) {
				fill(F[d], F[d] + (1 << n), Z(0));
				for (int i = 1; i < 1 << n; i <<= 1)
					for (int j = 0; j < 1 << n; j += i << 1)
						for (int S = j; S < j + i; S ++)
							F[d - 1][S | i] += F[d - 1][S];
				for (int S = 0; S < 1 << n; S ++)
					if (int c = __builtin_popcount(S); c <= d && d <= c * 2)
						for (int i = d; i; i --)
							F[d][S] += G[S][i] * F[d - i][S];
				for (int i = 1; i < 1 << n; i <<= 1)
					for (int j = 0; j < 1 << n; j += i << 1)
						for (int S = j; S < j + i; S ++)
							F[d][S | i] -= F[d][S];
				for (int S = 0; S < 1 << n; S ++)
					if (__builtin_popcount(S) == d) func(S, F[d][S]);
					else F[d][S] = 0;
			}
			for (int S = 0; S < 1 << n; S ++)
				f[S] = F[__builtin_popcount(S)][S];
			return;
		}
		template<class Func = void(*)(int, Z&)>
		static constexpr void Online_Convolution(const int &n, Z *f, const Z *g, const Z &ori, const Func &func = [](int, Z&) {}) {
			(n <= 11) ? Online_Convolution_Naive(n, f, g, ori, func)
								: Online_Convolution_FWT(n, f, g, ori, func);
			return;
		}
	
	public:
		constexpr SetPowerSeries(): vector<Z>() {}
		explicit constexpr SetPowerSeries(const int n): vector<Z>(n) {}
		explicit constexpr SetPowerSeries(const vector<Z> &a): vector<Z>(a) {}
		constexpr SetPowerSeries(const initializer_list<Z> &a): vector<Z>(a) {}
		template<class _InputIterator, class = _RequireInputIter<_InputIterator>>
		explicit constexpr SetPowerSeries(_InputIterator __first, _InputIterator __last): vector<Z>(__first, __last) {}
		template<class F = Z(*)(int)> explicit constexpr SetPowerSeries(int n, F f): vector<Z>(n) {
			for (int i = 0; i < n; i ++) (*this)[i] = f(i);
		}
		
		friend constexpr SetPowerSeries Subset_Zeta(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						f[S | i] += f[S];
			return f;
		}
		friend constexpr SetPowerSeries Subset_Mobius(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						f[S | i] -= f[S];
			return f;
		}
		friend constexpr SetPowerSeries Supset_Zeta(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						f[S] += f[S | i];
			return f;
		}
		friend constexpr SetPowerSeries Supset_Mobius(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++)
						f[S] -= f[S | i];
			return f;
		}
		friend constexpr SetPowerSeries FWT(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++) {
						Z x = f[S], y = f[S | i];
						f[S] = x + y, f[S | i] = x - y;
					}
			return f;
		}
		friend constexpr SetPowerSeries IFWT(SetPowerSeries f) {
			constexpr Z inv2 = (MOD + 1) >> 1;
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			for (int i = 1; i < 1 << n; i <<= 1)
				for (int j = 0; j < 1 << n; j += i << 1)
					for (int S = j; S < j + i; S ++) {
						Z x = f[S], y = f[S | i];
						f[S] = inv2 * (x + y), f[S | i] = inv2 * (x - y);
					}
			return f;
		}
		
		friend constexpr SetPowerSeries And_Convolution(SetPowerSeries f, SetPowerSeries g) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(static_cast<int>(g.size()) == 1 << n);
			f = Supset_Zeta(f), g = Supset_Zeta(g);
			for (int S = 0; S < 1 << n; S ++)
				f[S] *= g[S];
			return Supset_Mobius(f);
		}
		friend constexpr SetPowerSeries Or_Convolution(SetPowerSeries f, SetPowerSeries g) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(static_cast<int>(g.size()) == 1 << n);
			f = Subset_Zeta(f), g = Subset_Zeta(g);
			for (int S = 0; S < 1 << n; S ++)
				f[S] *= g[S];
			return Subset_Mobius(f);
		}
		friend constexpr SetPowerSeries Xor_Convolution(SetPowerSeries f, SetPowerSeries g) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(static_cast<int>(g.size()) == 1 << n);
			f = FWT(f), g = FWT(g);
			for (int S = 0; S < 1 << n; S ++)
				f[S] *= g[S];
			return IFWT(f);
		}
		
		friend constexpr SetPowerSeries Subset_Convolution(const SetPowerSeries &f, const SetPowerSeries &g) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == (1 << n));
			assert(static_cast<int>(g.size()) == (1 << n));
			SetPowerSeries ret(1 << n);
			Convolution(n, f.data(), g.data(), ret.data());
			return ret;
		}
		friend constexpr SetPowerSeries Subset_Division(SetPowerSeries f, const SetPowerSeries &g) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == (1 << n));
			assert(static_cast<int>(g.size()) == (1 << n));
			assert(g[0].x == 1);
			Division(n, f.data(), g.data());
			return f;
		}
		template<class Func = void(*)(int, Z&)>
		friend constexpr SetPowerSeries Online_Subset_Convolution(const SetPowerSeries &g, const Z &ori, const Func &func = [](int, Z&) {}) {
			const int n = __lg(g.size());
			assert(0 <= n && n <= N);
			assert(static_cast<int>(g.size()) == 1 << n);
			SetPowerSeries f(1 << n);
			Online_Convolution(n, f.data(), g.data(), ori, func);
			return f;
		}
		
		friend constexpr SetPowerSeries Inv(const SetPowerSeries &f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(f[0].x == 1);
			SetPowerSeries ret(1 << n);
			ret[0] = 1, Division(n, ret.data(), f.data());
			return ret;
		}
		friend constexpr SetPowerSeries Exp(const SetPowerSeries &f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(f[0].x == 0);
			assert(static_cast<int>(f.size()) == (1 << n));
			SetPowerSeries g(1 << n); g[0] = 1;
			for (int k = 0; k < n; k ++)
				Convolution(k, g.data(), f.data() + (1 << k), g.data() + (1 << k));
			return g;
		}
		friend constexpr SetPowerSeries Log(SetPowerSeries f) {
			const int n = __lg(f.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(f[0].x == 1);
			SetPowerSeries g(1 << n);
			for (int k = n - 1; ~k; k --) {
				copy(f.begin() + (1 << k), f.begin() + (1 << (k + 1)), g.begin() + (1 << k));
				Division(k, g.data() + (1 << k), f.data());
			}
			return g;
		}
		
		template<class EGF>
		friend constexpr SetPowerSeries EGF_Composite(const EGF &egf, const SetPowerSeries &f) {
			const int n = __lg(f.size());
			assert(static_cast<int>(egf.size()) == n + 1);
			assert(static_cast<int>(f.size()) == 1 << n);
			assert(f[0].x == 0);
			SetPowerSeries g(1 << n);
			Composite(n, egf.data(), f.data(), g.data());
			return g;
		}
		template<class Poly>
		friend constexpr SetPowerSeries Poly_Composite(Poly p, SetPowerSeries f) {
			const int n = __lg(f.size()), m = static_cast<int>(p.size());
			assert(static_cast<int>(f.size()) == 1 << n);
			if (!m) return SetPowerSeries(1 << n);
			vector<Z> egf(n + 1);
			for (int i = 0; i <= n; i ++) {
				Z x = 0;
				for (int j = m; -- j >= 0; )
					(x *= f[0]) += p[j];
				egf[i] = x;
				for (int j = 1; j < m; j ++)
					p[j - 1] = p[j] * j;
				p[m - 1] = 0;
			}
			f[0] = 0;
			return EGF_Composite(egf, f);
		}
		
		constexpr SetPowerSeries& operator += (const SetPowerSeries &rhs) {
			const int n = __lg(rhs.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(this -> size()) == 1 << n);
			for (int S = 0; S < 1 << n; S ++) (*this)[S] += rhs[S];
			return *this;
		}
		constexpr SetPowerSeries& operator -= (const SetPowerSeries &rhs) {
			const int n = __lg(rhs.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(this -> size()) == 1 << n);
			for (int S = 0; S < 1 << n; S ++) (*this)[S] -= rhs[S];
			return *this;
		}
		constexpr SetPowerSeries& operator *= (const SetPowerSeries &rhs) {
			const int n = __lg(rhs.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(this -> size()) == 1 << n);
			return (*this) = Subset_Convolution(*this, rhs);
		}
		constexpr SetPowerSeries& operator /= (const SetPowerSeries &rhs) {
			const int n = __lg(rhs.size());
			assert(0 <= n), assert(n <= N);
			assert(static_cast<int>(this -> size()) == 1 << n);
			return (*this) = Subset_Division(*this, rhs);
		}
		
		friend constexpr SetPowerSeries operator + (SetPowerSeries f, const SetPowerSeries &g) {
			return f += g;
		}
		friend constexpr SetPowerSeries operator - (SetPowerSeries f, const SetPowerSeries &g) {
			return f -= g;
		}
		friend constexpr SetPowerSeries operator * (SetPowerSeries f, const SetPowerSeries &g) {
			return f *= g;
		}
		friend constexpr SetPowerSeries operator / (SetPowerSeries f, const SetPowerSeries &g) {
			return f /= g;
		}
}; using SPS = SetPowerSeries<>;
```