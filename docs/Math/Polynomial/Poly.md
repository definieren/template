```cpp
namespace Polynomial {
using u64 = unsigned long long;

constexpr unsigned Mod = Z::GetMod();
constexpr size_t logord = __builtin_ctz(Mod - 1);
constexpr Z g = 3, zeta = g.Pow((Mod - 1) >> logord), inv_zeta = zeta.inv();
vector<Z> _root{Z::raw(1)}, _inv_root{Z::raw(1)};

Z* root(const size_t &n) {
	const size_t sz = _root.size();
	if (sz < n) {
		_root.resize(n);
		for (size_t i = __builtin_ctz(sz); static_cast<size_t>(1) << i < n; i ++) {
			const size_t j = 1 << i;
			_root[j] = zeta.Pow(1 << (logord - i - 2));
			for (size_t k = j + 1; k < j * 2; k ++)
				_root[k] = _root[k - j] * _root[j];
		}
	}
	return _root.data();
}
Z* inv_root(const size_t &n) {
	const size_t sz = _inv_root.size();
	if (sz < n) {
		_inv_root.resize(n);
		for (size_t i = __builtin_ctz(sz); static_cast<size_t>(1) << i < n; i ++) {
			const size_t j = 1 << i;
			_inv_root[j] = inv_zeta.Pow(1 << (logord - i - 2));
			for (size_t k = j + 1; k < j * 2; k ++)
				_inv_root[k] = _inv_root[k - j] * _inv_root[j];
		}
	}
	return _inv_root.data();
}
void DFT(vector<Z> &a) {
	const size_t n = a.size();
	static vector<u64> b; b.resize(n), b[0] = 0;
	for (size_t j = 0, k = n >> 1; j < k; j ++) {
		u64 u = (u64)a[j], v = (u64)a[j | k];
		b[j] = u + v, b[j | k] = u + Mod - v;
	}
	Z* rt = root(n >> 1);
	for (size_t k = n >> 2; k >= 1; k >>= 1) {
		for (size_t j = 0; j < k; j ++) {
			u64 u = b[j], v = b[j | k] % Mod;
			b[j] = u + v, b[j | k] = u + Mod - v;
		}
		for (size_t i = k << 1, m = 1; i < n; i += k << 1, m ++)
			for (size_t j = 0; j < k; j ++) {
				u64 u = b[i | j], v = b[i | j | k] * (u64)rt[m] % Mod;
				b[i | j] = u + v, b[i | j | k] = u + Mod - v;	
			}
	}
	for (size_t i = 0; i < n; i ++) a[i] = b[i];
	return;
}
void IDFT(vector<Z> &a) {
	const size_t n = a.size();
	Z* rt = inv_root(n / 2);
	for (size_t k = 1; k << 1 < n; k <<= 1) {
		for (size_t j = 0; j < k; j ++) {
			Z u = a[j], v = a[j | k];
			a[j] = u + v, a[j | k] = u - v;
		}
		for (size_t i = k << 1, m = 1; i < n; i += k << 1, m ++)
			for (size_t j = 0; j < k; j ++) {
				Z u = a[i | j], v = a[i | j | k];
				a[i | j] = u + v, a[i | j | k] = (u - v) * rt[m];
			}
	}
	Z inv = Z::raw(Mod - Mod / n);
	for (size_t j = 0, k = n >> 1; j < k; j ++) {
		Z u = a[j] * inv, v = a[j | k] * inv;
		a[j] = u + v, a[j | k] = u - v;
	}
	return;
}

struct Poly : public vector<Z> {
	Poly(): vector<Z>() { return; }
	explicit Poly(int n): vector<Z>(n) { return; }
	explicit Poly(const vector<Z> &a): vector<Z>(a) { return; }
	explicit Poly(const initializer_list<Z> &a): vector<Z>(a) { return; }
	template<class _InputIterator, class = _RequireInputIter<_InputIterator>>
	explicit constexpr Poly(_InputIterator __first, _InputIterator __last): vector<Z>(__first, __last) { return; }
	template<class F> explicit constexpr Poly(int n, F f): vector<Z>(n) {
		for (int i = 0; i < n; i ++) (*this)[i] = f(i);
		return;
	}
	
	Poly Shift(const int &k) const {
		if (k >= 0) { auto ret = *this; ret.insert(ret.begin(), k, Z(0)); return ret; }
		else return (int)this -> size() <= - k ? Poly() : Poly(this -> begin() - k, this -> end());
	}
	Poly Trunc(const size_t &k) const {
		Poly f = *this; f.resize(k); return f;
	}
	
	friend Poly operator - (Poly rhs) {
		const size_t n = rhs.size();
		for (size_t i = 0; i < n; i ++) rhs[i] = - rhs[i];
		return rhs;
	}
	friend Poly operator + (const Poly &rhs) {
		return rhs;
	}
	Poly& operator += (const Poly &rhs) {
		const size_t n = size(), m = rhs.size();
		this -> resize(std::max(n, m));
		for (size_t i = 0; i < m; i ++) (*this)[i] += rhs[i];
		return *this;
	}
	Poly& operator += (const Z &rhs) {
		if (__builtin_expect(empty(), 0))
			return (*this) = Poly{rhs};
		return (*this)[0] += rhs, *this;
	}
	Poly& operator -= (const Poly &rhs) {
		const size_t n = size(), m = rhs.size();
		this -> resize(std::max(n, m));
		for (size_t i = 0; i < m; i ++) (*this)[i] -= rhs[i];
		return *this;
	}
	Poly& operator -= (const Z &rhs) {
		if (__builtin_expect(empty(), 0))
			return (*this) = Poly{-rhs};
		return (*this)[0] -= rhs, *this;
	}
	Poly& operator *= (Poly rhs) {
		if (__builtin_expect(!size() || !rhs.size(), 0))
			return (*this) = Poly();
		if (size() < rhs.size()) std::swap(*this, rhs);
		const size_t n1 = size(), n2 = rhs.size();
		const size_t m = n1 + n2 - 1;
		if (n2 < static_cast<size_t>(1) << 7) {
			Poly ret(m);
			for (size_t i = 0; i < n1; i ++)
				for (size_t j = 0; j < n2; j ++)
					ret[i + j] += (*this)[i] * rhs[j];
			return *this = move(ret);
		}
		size_t n = 1;
		while (n <= m) n <<= 1;
		resize(n), DFT(*this), rhs.resize(n), DFT(rhs);
		for (size_t i = 0; i < n; i ++) (*this)[i] *= rhs[i];
		IDFT(*this), resize(m);
		return *this;
	}
	Poly& operator *= (const Z &rhs) {
		const size_t n = size();
		for (size_t i = 0; i < n; i ++) (*this)[i] *= rhs;
		return *this;
	}
	Poly& operator /= (const Poly &rhs) {
		const size_t n = size(), m = rhs.size();
		if (__builtin_expect(m > n, 0)) return (*this) = Poly();
		auto RG = rhs;
		reverse(begin(), end()), resize(n + 1 - m);
		reverse(RG.begin(), RG.end()), RG.resize(n + 1 - m);
		(*this) *= Inv(RG, n + 1 - m), (*this) = Trunc(n + 1 - m);
		reverse(begin(), end());
		return *this;
	}
	Poly& operator /= (const Z &rhs) {
		return (*this) *= rhs.inv();
	}
	friend Poly operator + (Poly lhs, const Poly &rhs) { return lhs += rhs; }
	friend Poly operator + (Poly lhs, const Z &rhs) { return lhs += rhs; }
	friend Poly operator + (const Z &lhs, Poly rhs) { return rhs += lhs; }
	friend Poly operator - (Poly lhs, const Poly &rhs) { return lhs -= rhs; }
	friend Poly operator - (Poly lhs, const Z &rhs) { return lhs -= rhs; }
	friend Poly operator - (const Z &lhs, Poly rhs) { return -(rhs -= lhs); }
	friend Poly operator * (Poly lhs, const Poly &rhs) { return lhs *= rhs; }
	friend Poly operator * (Poly lhs, const Z &rhs) { return lhs *= rhs; }
	friend Poly operator * (const Z &lhs, Poly rhs) { return rhs *= lhs; }
	friend Poly operator / (Poly lhs, const Poly &rhs) { return lhs /= rhs; }
	friend Poly operator / (Poly lhs, const Z &rhs) { return lhs /= rhs; }
	friend Poly Deriv(const Poly &f) {
		const size_t n = f.size();
		if (n == 0) return f;
		Poly g(n - 1);
		for (size_t i = 0; i + 1 < n; i ++)
			g[i] = (i + 1) * f[i + 1];
		return g;
	}
	friend Poly Integr(const Poly &f) {
		const size_t n = f.size();
		Poly g(n + 1);
		for (size_t i = 0; i < n; i ++)
			g[i + 1] = f[i] * comb.inv(i + 1);
		return g;
	}
	friend Poly Inv(const Poly &f, const int &lim) {
		assert(f.size() && f[0]);
		Poly g{f[0].inv()}; int k = 1;
		while (k < lim) {
			g.resize(k <<= 1);
			auto F = f.Trunc(k), G = g;
			DFT(F), DFT(G);
			for (int i = 0; i < k; i ++) F[i] = F[i] * G[i];
			IDFT(F), -- F[0];
			for (int i = 0; i < k / 2; i ++) F[i] = 0;
			DFT(F);
			for (int i = 0; i < k; i ++) F[i] *= G[i];
			IDFT(F);
			for (int i = k / 2; i < k; i ++) g[i] -= F[i];
		}
		return g.Trunc(lim);
	}
	friend pair<Poly, Poly> Div(const Poly &f, const Poly &g) {
		const size_t m = g.size();
		assert(m != 0);
		auto Q = f / g, R = (f - g * Q).Trunc(m - 1);
		return {Q, R};
	}
	friend Poly Div(const Poly &h, const Poly &f, const size_t &lim) {
		if (__builtin_expect(lim == 0, 0)) return Poly();
		if (__builtin_expect(h.empty(), 0)) return Poly(lim);
		assert(f.size() && f[0]);
		size_t k = 1;
		Poly g{f[0].inv()}, q{g[0] * h[0]}, F, G, H, Q;
		while (k < lim) {
			q.resize(k <<= 1);
			DFT(F = move(f.Trunc(k))), DFT(Q = q), H = move(h.Trunc(k));
			for (size_t i = 0; i < k; i ++) F[i] *= Q[i];
			IDFT(F);
			for (size_t i = 0; i < k / 2; i ++) F[i] = Z::raw(0);
			for (size_t i = k / 2; i < k; i ++) F[i] -= H[i];
			DFT(F), g.resize(k), DFT(G = g);
			for (size_t i = 0; i < k; i ++) F[i] *= G[i];
			IDFT(F);
			for (size_t i = k / 2; i < k; i ++) q[i] -= F[i];
			if (k < lim) {
				g.resize(k);
				DFT(F = f.Trunc(k));
				for (size_t i = 0; i < k; i ++) F[i] = F[i] * G[i];
				IDFT(F), -- F[0];
				for (size_t i = 0; i < k / 2; i ++) F[i] = Z::raw(0);
				DFT(F);
				for (size_t i = 0; i < k; i ++) F[i] *= G[i];
				IDFT(F);
				for (size_t i = k / 2; i < k; i ++) g[i] -= F[i];
			}
		}
		return q.Trunc(lim);
	}
	friend Poly Log(const Poly &f, const size_t &lim) {
		return Integr(Div(Deriv(f), f, lim)).Trunc(lim);
	}
	friend Poly Sqrt(const Poly &f, const size_t &lim) {
		if (__builtin_expect(f.empty(), 0)) return Poly(lim);
		assert(f[0] == 1); size_t k = 1;
		Poly g{1}, h{1}, F{1}, G, H;
		while (k < lim) {
			g.resize(k <<= 1);
			G = F, F = f.Trunc(k), H = h;
			for (size_t i = 0; i < k / 2; i ++) G[i] *= G[i];
			IDFT(G);
			for (size_t i = 0; i < k / 2; i ++) G[i] -= F[i] + F[i + k / 2];
			G.resize(k), H.resize(k), DFT(G), DFT(H);
			for (size_t i = 0; i < k; i ++) G[i] *= H[i];
			IDFT(G);
			for (size_t i = 0; i < k / 2; i ++) g[i + k / 2] -= G[i].div_2();
			if (k < lim) {
				h.resize(k), DFT(G = g), F = G;
				for (size_t i = 0; i < k; i ++) G[i] = G[i] * H[i];
				IDFT(G), -- G[0];
				for (size_t i = 0; i < k / 2; i ++) G[i] = Z::raw(0);
				DFT(G);
				for (size_t i = 0; i < k; i ++) G[i] *= H[i];
				IDFT(G);
				for (size_t i = k / 2; i < k; i ++) h[i] -= G[i];
			}
		}
		return g.Trunc(lim);
	}
	friend Poly Exp(const Poly &f, const size_t &lim) {
		if (__builtin_expect(f.empty(), 0)) {
			if (__builtin_expect(lim == 0, 0)) return Poly();
			Poly ret(lim); ret[0] = Z::raw(1); return ret;
		}
		assert(f[0] == 0);
		size_t k = 1;
		Poly g{1}, h{1}, F{1}, dF, G, dG, H;
		while (k < lim) {
			G = F, dG = move(Deriv(g).Trunc(k));
			dF = move(Deriv(f.Trunc(k)).Trunc(k)), DFT(F = dF);
			for (size_t i = 0; i < k; i ++) G[i] *= F[i];
			IDFT(G), g.resize(k <<= 1);
			for (size_t i = 0; i < k / 2; i ++) G[i] = dG[i] - G[i];
			G.insert(G.begin(), G.back()), G.back() = Z::raw(0);
			G.resize(k), DFT(G), h.resize(k), H = h, DFT(H);
			for (size_t i = 0; i < k; i ++) G[i] *= H[i];
			IDFT(G), G = G.Shift(k / 2 - 1), G.resize(k - 1);
			G = Integr(G + dF) - f.Trunc(k);
			DFT(G), H = move(g.Trunc(k)), DFT(H);
			for (size_t i = 0; i < k; i ++) G[i] *= H[i];
			IDFT(G);
			for (size_t i = k / 2; i < k; i ++) g[i] -= G[i];
			if (k < lim) {
				h.resize(k);
				DFT(G = g), DFT(H = h), F = G;
				for (size_t i = 0; i < k; i ++) G[i] = G[i] * H[i];
				IDFT(G), -- G[0];
				for (size_t i = 0; i < k / 2; i ++) G[i] = Z::raw(0);
				DFT(G);
				for (size_t i = 0; i < k; i ++) G[i] *= H[i];
				IDFT(G);
				for (size_t i = k / 2; i < k; i ++) h[i] -= G[i];
			}
		}
		return g.Trunc(lim);
	}
	friend Poly Pow(const Poly &f, const u64 &k, const size_t &lim) {
		size_t i = 0;
		while (i < f.size() && !f[i]) i ++;
		if (i == f.size() || i * k >= lim) return Poly(lim);
		Z v = f[i]; Poly g = f.Shift(- i) * v.inv();
		return Exp(Log(move(g), lim - i * k) * k, lim - i * k).Shift(i * k) * v.Pow(k);
	}
	friend Poly Pow(const Poly &f, const u64 &k1, const u64 &k2, const size_t &lim) {
		size_t i = 0;
		while (i < f.size() && !f[i]) i ++;
		if (i == f.size() || i * k1 >= lim) return Poly(lim);
		Z v = f[i]; Poly g = f.Shift(- i) * v.inv();
		return Exp(Log(move(g), lim - i * k1) * k1, lim - i * k1).Shift(i * k1) * v.Pow(k2);
	}
};
}
using Polynomial::DFT;
using Polynomial::IDFT;
using Polynomial::Poly;

Poly Berlekamp_Massey(const Poly &a) {
	Poly c, oldc; int f = - 1;
	for (int i = 0; i < (int)a.size(); i ++) {
		Z delta = a[i];
		for (int j = 0; j < (int)c.size(); j ++)
			delta -= c[j] * a[i - j - 1];
		if (!delta) continue;
		if (!~f) { c.resize(i + 1), f = i; continue; }
		Poly d = oldc; Z df = 0;
		d *= - 1, d.insert(d.begin(), 1);
		for (int j = 0; j < (int)d.size(); j ++)
			df += d[j] * a[f - j];
		(d *= delta / df).insert(d.begin(), i - f - 1, Z(0));
		auto tmp = c; c += d;
		if (i - tmp.size() > f - oldc.size())
			oldc = tmp, f = i;
	}
	c *= - 1, c.insert(c.begin(), 1);
	return c;
}
Z Linear_Recurrence(Poly F, Poly G, i64 n) {
	Poly P, Q;
	const int m = G.size();
	while (n) {
		Q = G;
		for (int i = 1; i < m; i += 2) Q[i] = - Q[i];
		P = F * Q, Q = G * Q;
		for (int i = 0; i < m; i ++) G[i] = Q[i << 1];
		for (int i = 0; i + 1 < m; i ++) F[i] = P[i << 1 | (n & 1)];
		n >>= 1;
	}
	return F[0] / G[0];
}
```