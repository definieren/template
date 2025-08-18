暂未完工。

还剩多点求值、快速插值、多项式复合、线性递推没写。

```cpp
template<class Z> 
struct polynomial: public std::vector<Z> {
private:
	static constexpr Z g = 3;
	static constexpr auto Mod = Z::getMod();
	static constexpr int log_ord = []() {
		auto x = Mod - 1;
		int y = 0;
		while (~x & 1) {
			x >>= 1, ++ y;
		}
		return y;
	}();
	static constexpr std::array<Z, log_ord + 1> invn = []() {
		std::array<Z, log_ord + 1> inv{};
		for (int i = 0; i <= log_ord; i ++) {
			inv[i] = Z(1 << i).inv();
		}
		return inv;
	}();
	static std::pair<Z*, Z*> get_root(const int &n) {
		static std::vector<Z> root{Z::from_raw(1)};
		static std::vector<Z> inv_root{Z::from_raw(1)};
		if (static_cast<int>(root.size()) < n) {
			int i = root.size();
			root.resize(n), inv_root.resize(n);
			for (; i != n; i <<= 1) {
				const Z w = g.pow(Mod / (i << 2)), iw = w.inv();
				for (int j = 0; j != i; j ++) {
					root[i + j] = root[j] * w;
					inv_root[i + j] = inv_root[j] * iw;
				}
			}
		}
		return {root.data(), inv_root.data()};
	}
	static constexpr int get_len(int n) {
		return n < 3 ? n : 2 << std::__lg(n - 1);
	}
	static void dif_n(Z *f, const int &n) {
		const Z* rt = get_root(n).first;
		for (int i = n; i >>= 1; ) {
			for (int j = 0, k = 0; j != n; j += i << 1, ++ k) {
				for (int p = j, q = j + i; p != j + i; ++ p, ++ q) {
					const Z u = f[p], v = f[q] * rt[k];
					f[p] = u + v, f[q] = u - v;
				}
			}
		}
	}
	static void dit_n(Z *f, const int &n) {
		const Z* irt = get_root(n).second;
		for (int i = 1; i != n; i <<= 1) {
			for (int j = 0, k = 0; j != n; j += i << 1, ++ k) {
				for (int p = j, q = j + i; p != j + i; ++ p, ++ q) {
					const Z u = f[p], v = f[q];
					f[p] = u + v, f[q] = (u - v) * irt[k];
				}
			}
		}
		const Z inv = invn[std::__lg(n)];
		for (int i = 0; i < n; i ++) {
			f[i] *= inv;
		}
	}
	static void dif_rhalf_n(Z *f, const int &n) {
		const Z* rt = get_root(n).first;
		for (int i = n, m = 1; i >>= 1; m <<= 1) {
			for (int j = 0, k = 0; j != n; j += i << 1, ++ k) {
				for (int p = j, q = j + i; p != j + i; ++ p, ++ q) {
					const Z u = f[p], v = f[q] * rt[m + k];
					f[p] = u + v, f[q] = u - v; 
				}
			}
		}
	}
	static void dit_rhalf_n(Z *f, const int &n) {
		const Z* irt = get_root(n).second;
		for (int i = 1, m = n; m >>= 1; i <<= 1) {
			for (int j = 0, k = 0; j != n; j += i << 1, ++ k) {
				for (int p = j, q = j + i; p != j + i; ++ p, ++ q) {
					const Z u = f[p], v = f[q];
					f[p] = u + v, f[q] = (u - v) * irt[m + k];
				}
			}
		}
		const Z inv = invn[std::__lg(n)];
		for (int i = 0; i < n; i ++) {
			f[i] *= inv;
		}
	}
	static void neg_n(const Z *a, const int &n, Z *b) {
		for (int i = 0; i < n; i ++) {
			b[i] = -a[i];
		}
	}
	static void add_n(const Z *a, const Z *b, const int &n, Z *c) {
		for (int i = 0; i < n; i ++) {
			c[i] = a[i] + b[i];
		}
	}
	static void sub_n(const Z *a, const Z *b, const int &n, Z *c) {
		for (int i = 0; i < n; i ++) {
			c[i] = a[i] - b[i];
		}
	}
	static void dot_n(const Z *a, const Z *b, const int &n, Z *c) {
		for (int i = 0; i < n; i ++) {
			c[i] = a[i] * b[i];
		}
	}
	static void mul_c_n(const Z *a, const Z &c, const int &n, Z *b) {
		for (int i = 0; i < n; i ++) {
			b[i] = a[i] * c;
		}
	}
	static void deriv_n(const Z *a, const int &n, Z *b) {
		for (int i = 1; i != n; i ++) {
			b[i - 1] = a[i] * i;
		}
		b[n - 1] = Z::from_raw(0);
	}
	static void integ_n(const Z *a, const int &n, Z *b) {
		comb.init(n);
		for (int i = n - 1; i; i --) {
			b[i] = a[i - 1] * comb.inv(i);
		}
		b[0] = Z::from_raw(0);
	}
	static void zero_n(Z *a, const int &n) {
		for (int i = 0; i < n; i ++) {
			a[i] = Z::from_raw(0);
		}
	}
	static void copy_n(const Z *a, const int &n, Z *b) {
		memcpy(b, a, sizeof(Z) * n);
	}
	static polynomial convolution_fft(polynomial f, polynomial g) {
		const int n = f.size() + g.size() - 1, m = get_len(n);
		f.resize(m), dif_n(f.data(), m);
		g.resize(m), dif_n(g.data(), m);
		dot_n(f.data(), g.data(), m, f.data());
		dit_n(f.data(), m), f.resize(n);
		return f;
	}
	static polynomial convolution_naive(const polynomial &f, const polynomial &g) {
		const int n = f.size(), m = g.size();
		if (__builtin_expect(f.empty() || g.empty(), 0)) {
			return polynomial{};
		}
		polynomial fg(n + m - 1);
		for (int i = 0; i != n; i ++) {
			for (int j = 0; j != m; j ++) {
				fg[i + j] += f[i] * g[j];
			}
		}
		return fg;
	}
	polynomial pow_one(const int &n, const int &k) const {
		if (__builtin_expect(k == 0, 0)) {
			polynomial f(this->size());
			f[0] = 1;
			return f;
		}
		return (k == 1) ? mod_xk(n) : (mod_xk(n).log() * k).exp();
	}
	polynomial pow_ord_zero(const int &n, int k_mod_p, int k_mod_phi_p) const {
		if ((*this)[0] == Z::from_raw(1)) {
			return pow_one(n, k_mod_p);
		}
		return (*this * (*this)[0].inv()).pow_one(n, k_mod_p) * (*this)[0].pow(k_mod_phi_p);
	}
	polynomial sqrt_ord_zero(const Z &s) const {
		constexpr Z c = -Z(2).inv();
		const int shrink_len = this->size();
		const int n = get_len(shrink_len);
		polynomial res(shrink_len), inv_res(shrink_len);
		polynomial dft_res(n), dft_inv_res(n), f(n);
		int N = 1, N2 = 2;
		res[0] = dft_res[0] = s, inv_res[0] = s.inv();
		while (N < shrink_len) {
			const int newN = (N2 == n) ? shrink_len : N2;
			dot_n(dft_res.data(), dft_res.data(), N, dft_res.data()), dit_n(dft_res.data(), N);
			sub_n(dft_res.data(), this->data(), N, dft_res.data() + N);
			sub_n(dft_res.data() + N, this->data() + N, newN - N, dft_res.data() + N);
			zero_n(dft_res.data(), N), dif_n(dft_res.data(), N2);
			copy_n(inv_res.data(), N, dft_inv_res.data()), dif_n(dft_inv_res.data(), N2);
			dot_n(dft_res.data(), dft_inv_res.data(), N2, dft_res.data()), dit_n(dft_res.data(), N2);
			mul_c_n(dft_res.data() + N, c, newN - N, res.data() + N);
			if (__builtin_expect(N2 < n, 1)) {
				copy_n(res.data(), N2, f.data()), dif_n(f.data(), N2);
				copy_n(f.data(), N2, dft_res.data());
				dot_n(f.data(), dft_inv_res.data(), N2, f.data()), dit_n(f.data(), N2);
				zero_n(f.data(), N), dif_n(f.data(), N2);
				dot_n(f.data(), dft_inv_res.data(), N2, f.data()), dit_n(f.data(), N2);
				neg_n(f.data() + N, N, inv_res.data() + N);
			}
			N <<= 1, N2 <<= 1;
		}
		return res;
	}
public:
	polynomial(): std::vector<Z>() {}
	explicit polynomial(const int &n): std::vector<Z>(n) {}
	explicit polynomial(const std::vector<Z> &a): std::vector<Z>(a) {}
	explicit polynomial(const std::initializer_list<Z> &a): std::vector<Z>(a) {}
	template<class _InputIterator, class = std::_RequireInputIter<_InputIterator>>
	explicit polynomial(_InputIterator __first, _InputIterator __last): std::vector<Z>(__first, __last) {}
	template<class F, class = Z(*)(int)>
	explicit polynomial(const int &n, const F &f): std::vector<Z>(n) {
		for (int i = 0; i != n; i ++) {
			(*this)[i] = f(i);
		}
	}
	
	int ord() const {
		int ord = 0;
		while (ord != static_cast<int>(this->size()) && !(*this)[ord]) {
			++ ord;
		}
		return ord;
	}
	int deg() const {
		int deg = static_cast<int>(this->size()) - 1;
		while (deg >= 0 && !(*this)[deg]) {
			-- deg;
		}
		return deg;
	}
	void remove0() {
		while (this->size() && !this->back()) {
			this->pop_back();
		}
	}
	polynomial mul_xk(const int &k) const {
		auto f = *this;
		f.insert(f.begin(), k, Z::from_raw(0));
		return f;
	}
	polynomial div_xk(const int &k) const {
		return polynomial(this->begin() + std::min(static_cast<int>(this->size()), k), this->end());
	}
	polynomial mod_xk(const int &k) const {
		return polynomial(this->begin(), this->begin() + std::min(static_cast<int>(this->size()), k));
	}
	polynomial compose_cx(const Z &c) const {
		const int n = this->size();
		auto f = *this;
		Z ci = Z::from_raw(1);
		for (int i = 1; i < n; i ++) {
			f[i] *= (ci *= c);
		}
		return f;
	}
	Z operator () (const Z &x) const {
		const int n = this->size();
		Z y = Z::from_raw(0);
		for (int i = n - 1; i >= 0; i --) {
			(y *= x) += (*this)[i];
		}
		return y;
	}
	polynomial operator + () const {
		return *this;
	}
	polynomial operator - () const {
		auto f = *this;
		neg_n(f.data(), f.size(), f.data());
		return f;
	}
	polynomial& operator += (const Z &rhs) {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial(1, rhs);
		}
		return (*this)[0] += rhs, *this;
	}
	polynomial& operator -= (const Z &rhs) {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial(1, -rhs);
		}
		return (*this)[0] -= rhs, *this;
	}
	polynomial& operator *= (const Z &rhs) {
		mul_c_n(this->data(), rhs, this->size(), this->data());
		return *this;
	}
	polynomial& operator /= (const Z &rhs) {
		mul_c_n(this->data(), rhs.inv(), this->size(), this->data());
		return *this;
	}
	polynomial& operator += (const polynomial &rhs) {
		if (this->size() < rhs.size()) {
			this->resize(rhs.size());
		}
		add_n(this->data(), rhs.data(), rhs.size(), this->data());
		return *this;
	}
	polynomial& operator -= (const polynomial &rhs) {
		if (this->size() < rhs.size()) {
			this->resize(rhs.size());
		}
		sub_n(this->data(), rhs.data(), rhs.size(), this->data());
		return *this;
	}
	polynomial& operator *= (const polynomial &rhs) {
		return *this = *this * rhs;
	}
	friend polynomial operator + (polynomial lhs, const Z &rhs) {
		return lhs += rhs;
	}
	friend polynomial operator - (polynomial lhs, const Z &rhs) {
		return lhs -= rhs;
	}
	friend polynomial operator * (polynomial lhs, const Z &rhs) {
		return lhs *= rhs;
	}
	friend polynomial operator / (polynomial lhs, const Z &rhs) {
		return lhs /= rhs;
	}
	friend polynomial operator + (polynomial lhs, const polynomial &rhs) {
		return lhs += rhs;
	}
	friend polynomial operator - (polynomial lhs, const polynomial &rhs) {
		return lhs -= rhs;
	}
	friend polynomial operator * (const polynomial &f, const polynomial &g) {
		if (std::min(f.size(), g.size()) > 8 && std::max(f.size(), g.size()) > 128) {
			return convolution_fft(f, g);
		} else {
			return convolution_naive(f, g);
		}
	}
	polynomial deriv() const {
		if (__builtin_expect(this->emtpy())) {
			return *this;
		}
		auto f = *this;
		deriv_n(f.data(), f.size(), f.data());
		return f.pop_back(), f;
	}
	polynomial integ() const {
		auto f = *this;
		f.resize(this->size() + 1);
		integ_n(f.data(), f.size(), f.data());
		return f.pop_back(), f;
	}
	polynomial inv() const {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial{};
		}
		assert((*this)[0] != Z::from_raw(0));
		const int shrink_len = this->size();
		const int n = get_len(shrink_len);
		polynomial res(shrink_len), f(n), g(n);
		res[0] = (*this)[0].inv();
		int N = 1, N2 = 2;
		while (N < shrink_len) {
			const int newN = (N2 == n) ? shrink_len : N2;
			copy_n(this->data(), newN, f.data()), dif_n(f.data(), N2);
			copy_n(res.data(), N, g.data()), dif_n(g.data(), N2);
			dot_n(f.data(), g.data(), N2, f.data()), dit_n(f.data(), N2);
			zero_n(f.data(), N), dif_n(f.data(), N2);
			dot_n(f.data(), g.data(), N2, f.data()), dit_n(f.data(), N2);
			neg_n(f.data() + N, newN - N, res.data() + N);
			N = newN, N2 = N << 1;
		}
		return res;
	}
	polynomial div(const polynomial &g) const {
		if (g.size() > this->size()) {
			return polynomial{};
		}
		assert(g.size());
		const int n = this->size() - g.size() + 1;
		polynomial q(n), r(n);
		std::reverse_copy(this->begin() + g.size() - 1, this->end(), q.begin());
		std::reverse_copy(g.begin() + std::max(0, static_cast<int>(g.size()) - n), g.end(), r.begin());
		q = (q * r.inv()).mod_xk(n), std::reverse(q.begin(), q.end());
		return q;
	}
	std::pair<polynomial, polynomial> div_mod(const polynomial &g) const {
		if (g.size() > this->size()) {
			return {polynomial{}, *this};
		}
		assert(g.size());
		const int n = g.size() - 1;
		auto q = div(g);
		return {q, mod_xk(n) - (g.mod_xk(n) * q.mod_xk(n)).mod_xk(n)};
	}
	polynomial exp() const {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial{};
		}
		assert((*this)[0] == Z::from_raw(0));
		if (__builtin_expect(this->size() == 1, 0)) {
			return polynomial{Z::from_raw(1)};
		}
		const int shrink_len = this->size();
		const int n = get_len(shrink_len);
		int N = 1, N2 = 2;
		polynomial res(shrink_len), inv_res(shrink_len);
		polynomial dft_res(n), dft_inv_res(n), f(n);
		res[0] = inv_res[0] = dft_res[0] = dft_res[1] =1;
		while (N < shrink_len) {
			const int newN = (N2 == n) ? shrink_len : N2;
			deriv_n(this->data(), N, f.data()), dif_n(f.data(), N);
			dot_n(f.data(), dft_res.data(), N, f.data()), dit_n(f.data(), N);
			f[N - 1] = -f[N - 1], deriv_n(res.data(), N, f.data() + N);
			sub_n(f.data() + N, f.data(), N - 1, f.data() + N);
			zero_n(f.data(), N - 1), dif_n(f.data(), N2);
			copy_n(inv_res.data(), N, dft_inv_res.data()), dif_n(dft_inv_res.data(), N2);
			dot_n(f.data(), dft_inv_res.data(), N2, f.data()), dit_n(f.data(), N2);
			integ_n(f.data(), N2, f.data()), zero_n(f.data(), N);
			sub_n(f.data() + N, this->data() + N, newN - N, f.data() + N), dif_n(f.data(), N2);
			dot_n(f.data(), dft_res.data(), N2, f.data()), dit_n(f.data(), N2);
			neg_n(f.data() + N, newN - N, res.data() + N);
			if (__builtin_expect(N2 < n, 1)) {
				copy_n(res.data(), N2, dft_res.data()), dif_n(dft_res.data(), N2 + N2);
				copy_n(dft_res.data(), N2, f.data());
				dot_n(f.data(), dft_inv_res.data(), N2, f.data()), dit_n(f.data(), N2);
				zero_n(f.data(), N), dif_n(f.data(), N2);
				dot_n(f.data(), dft_inv_res.data(), N2, f.data()), dit_n(f.data(), N2);
				neg_n(f.data() + N, N, inv_res.data() + N);
			}
			N <<= 1, N2 <<= 1;
		}
		return res;
	}
	polynomial log() const {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial{};
		}
		assert((*this)[0] == Z::from_raw(1));
		polynomial f(this->size());
		deriv_n(this->data(), this->size(), f.data());
		f *= inv(), f.resize(this->size());
		integ_n(f.data(), this->size(), f.data());
		return f;
	}
	polynomial pow(int k_chkmin_sz, int k_mod_p, int k_mod_phi_p) const {
		if (__builtin_expect(this->empty(), 0)) {
			return polynomial{};
		}
		if (__builtin_expect(k_chkmin_sz == 0, 0)) {
			polynomial f(this->size());
			return f[0] = Z::from_raw(1), f;
		}
		const int i = ord();
		if (static_cast<i64>(i) * k_chkmin_sz >= static_cast<i64>(this->size())) {
			return polynomial(this->size());
		}
		return div_xk(i).pow_ord_zero(this->size() - i * k_chkmin_sz, k_mod_p, k_mod_phi_p).mul_xk(i * k_chkmin_sz);
	}
	template<class T> inline polynomial pow(T x) const {
		if (x < 0) return inv().pow(-x);
		return pow(x < static_cast<T>(this->size()) ? x : this->size(), x % Mod, x % (Mod - 1));
	}
	std::pair<bool, polynomial> sqrt() const {
		const int i = ord();
		if (i == static_cast<int>(this->size())) {
			return {true, polynomial(this->size())};
		}
		if (i & 1) {
			return {false, polynomial{}};
		}
		auto [o, s] = (*this)[i].sqrt();
		if (o == false) {
			return {false, polynomial{}};
		}
		auto x = div_xk(i);
		x.insert(x.end(), i >> 1, 0);
		return {true, x.sqrt_ord_zero(s).mul_xk(i >> 1)};
	}
	polynomial composional_inv() const {
		assert(this->size() >= 2U);
		assert((*this)[0] == Z::from_raw(0));
		assert((*this)[1] != Z::from_raw(0));
		constexpr Z inv2 = Z(2).inv();
		const int shrink_len = this->size();
		const int n = get_len(shrink_len << 1);
		const Z *irt = get_root(2 * n).second;
		const Z invf1 = (*this)[1].inv();
		
		polynomial P(n), Q(n);
		polynomial dftP(n << 1), dftQ(n << 1);
		P[0] = Z::from_raw(1);
		Z v = -1;
		for (int i = 1; i < shrink_len; i ++) {
			v *= invf1, Q[i] = (*this)[i] * v;
		}
		int N = shrink_len, M = 1, newN = 0;
		for (; N > 1; N = newN, M <<= 1) {
			newN = (N + 1) >> 1;
			const int len = get_len((N * M) << 2);
			zero_n(dftP.data(), len), zero_n(dftQ.data(), len);
			for (int j = 0; j < M; j ++) {
				copy_n(P.data() + j * N, N, dftP.data() + 2 * j * N);
				copy_n(Q.data() + j * N, N, dftQ.data() + 2 * j * N);
			}
			
			dftQ[2 * N * M] = 1;
			dif_n(dftP.data(), len);
			dif_n(dftQ.data(), len);
			P.resize(len >> 1), Q.resize(len >> 1);
			for (int i = 0; i < len; i += 2) {
				Q[i >> 1] = dftQ[i] * dftQ[i + 1];
			}
			if (N & 1) {
				for (int i = 0; i < len; i += 2) {
					P[i >> 1] = (dftP[i] * dftQ[i + 1] + dftP[i + 1] * dftQ[i]) * inv2;
				}
			} else {
				for (int i = 0; i < len; i += 2) {
					P[i >> 1] = (dftP[i] * dftQ[i + 1] - dftP[i + 1] * dftQ[i]) * irt[i >> 1] * inv2;
				}
			}
			dit_n(P.data(), len >> 1);
			dit_n(Q.data(), len >> 1);
			if (N * M * 4 == len) {
				-- Q[0];
			}
			
			for (int j = 1; j < M * 2; j ++) {
				copy_n(P.data() + j * N, newN, P.data() + j * newN);
				copy_n(Q.data() + j * N, newN, Q.data() + j * newN);
			}
		}
		
		P.resize(M * newN);
		std::reverse(P.begin(), P.end());
		P.resize(shrink_len);
		for (int i = 1; i < shrink_len; i ++) {
			P[i] *= (shrink_len - 1) * comb.inv(i);
		}
		std::reverse(P.begin(), P.end());
		P = P.pow_one(shrink_len - 1, (int)(-comb.inv(shrink_len - 1))) * invf1;
		P.insert(P.begin(), 0);
		return P;
	}
};
using poly = polynomial<Z>;
```