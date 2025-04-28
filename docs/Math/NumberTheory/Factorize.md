# 质因数分解

`Is_Prime` 是 Miller-Rabin 判断是否是素数，`Factorize` 是 Pollard-Rho 质因数分解，返回一个 `vector<pair<u64, u32>>`，`first` 是质因数，`second` 是质因数幂次。

```cpp
template<class U0, class U1>
struct Montgomery {
	constexpr static unsigned B0 = sizeof(U0) * 8U;
	U0 n, nr, rs, np;
	
	constexpr Montgomery(const U0 &Mod) {
		SetMod(Mod);
	}
	
	constexpr U0 GetMod() const noexcept {
		return n;
	}
	constexpr void SetMod(const U0& Mod) {
		assert(Mod >= 2), assert(Mod % 2 == 1);
		assert((Mod >> (B0 - 2)) == 0);
		n = nr = Mod, rs = -static_cast<U1>(n) % n;
		for (u32 i = 0; i < __lg(B0); i ++) {
			nr *= 2 - n * nr;
		}
		np = Reduce(static_cast<U0>(1), rs);
	}
	constexpr U0 Reduce(const U0& x) const noexcept {
		const U0 q = x * nr;
		const U0 m = (static_cast<U1>(q) * n) >> B0;
		return n - m;
	}
	constexpr U0 Reduce(const U0& x, const U0& y) const noexcept {
		const U1 t = static_cast<U1>(x) * y;
		const U0 c = t, d = t >> B0;
		const U0 q = c * nr;
		const U0 m = (static_cast<U1>(q) * n) >> B0;
		return d + n - m;
	}
	constexpr U0 Reduce(const U0& x, const U0& y, const U0& z) const noexcept {
		const U1 t = static_cast<U1>(x) * y;
		const U0 c = t, d = t >> B0;
		const U0 q = c * nr;
		const U0 m = (static_cast<U1>(q) * n) >> B0;
		return z + d + n - m;
	}
	constexpr U0 val(const U0& x) const noexcept {
		const u64 t = Reduce(x);
		return (t == n) ? static_cast<U0>(0) : t;
	}
	constexpr U0 zero() const noexcept {
		return static_cast<U0>(0);
	}
	constexpr U0 one() const noexcept {
		return np;
	}
	constexpr U0 raw(const U0& x) const noexcept {
		return Reduce(x, rs);
	}
	template<class U> requires std::unsigned_integral<U>
	constexpr U0 trans(const U& x) const noexcept {
		if (__builtin_expect(x < n, 1)) {
			return raw(x);
		}
		return Reduce(x % n, rs);
	}
	template<class S> requires std::signed_integral<S>
	constexpr U0 trans(S x) const noexcept {
		if (__builtin_expect(0 <= x && x < static_cast<S>(n), 1)) {
			return Raw(x);
		}
		if ((x %= static_cast<S>(n)) < 0) {
			(x += static_cast<S>(n)) %= static_cast<S>(n);
		}
		return Reduce(x, rs);
	}
	constexpr U0 neg(const U0& x) const noexcept {
		return (x != 0) ? (2 * n - x) : x;
	}
	constexpr U0 inc(const U0& x) const noexcept {
		return add(x, np);
	}
	constexpr U0 dec(const U0& x) const noexcept {
		return sub(x, np);
	}
	constexpr U0 add(const U0& x, const U0& y) const noexcept {
		return (x + y >= 2 * n) ? (x + y - 2 * n) : (x + y);
	}
	constexpr U0 sub(const U0& x, const U0& y) const noexcept {
		return (x < y) ? (x - y + 2 * n) : (x - y);
	}
	constexpr U0 mul(const U0& x, const U0& y) const noexcept {
		return Reduce(x, y);
	}
	constexpr U0 mul_add(const U0& x, const U0& y, const U0& z) const noexcept {
		return Reduce(x, y, z);
	}
	constexpr bool same(const U0& x, const U0& y) const noexcept {
		const U0 dif = x - y;
		return (dif == 0) || (dif == n) || (dif == -n);
	}
};

constexpr bool Is_Prime(u64 x) noexcept {
	if (x <= 1) {
		return false;
	}
	if (x % 2 == 0) {
		return x == 2;
	}
	
	constexpr array<u64, 10> Base{2, 7, 61, 2, 325, 9375, 28178, 450775, 9780504, 1795265022};
	const u32 s = __builtin_ctzll(x - 1);
	const u64 d = (x - 1) >> s;
	const int q = 63 ^ __builtin_clzll(d);
	const Montgomery<u64, u128> Mod(x);
	const u32 l = (x >> 32) ? 3 : 0;
	const u32 r = (x >> 32) ? 10 : 3;
	for (u32 _ = l; _ < r; _ ++) {
		u64 base = Base[_];
		if (base % x == 0) {
			continue;
		}
		base = Mod.trans(base);
		u64 a = base;
		for (int i = q - 1; ~i; i --) {
			a = Mod.mul(a, a);
			if ((d >> i) & 1) {
				a = Mod.mul(a, base);
			}
		}
		if (Mod.same(a, Mod.one())) {
			continue;
		}
		for (u32 t = 1; t < s && !Mod.same(a, x - Mod.one()); ++ t) {
			a = Mod.mul(a, a);
		}
		if (!Mod.same(a, x - Mod.one())) {
			return false;
		}
	}
	return true;
}

u64 Find_Prime_Factor(const u64& n) {
	if (n % 2 == 0) {
		return 2ULL;
	}
	const Montgomery<u64, u128> Mod(n);
	const u64 C1 = 1, C2 = 2, M = 512;
	u64 Z1 = 1, Z2 = 2, ans = 0;
	auto find = [&]() {
		u64 z1 = Z1, z2 = Z2;
		for (u64 k = M; ; k *= 2) {
			const u64 x1 = z1 + n, x2 = z2 + n;
			for (u64 j = 0; j < k; j += M) {
				const u64 y1 = z1, y2 = z2;
				u64 q1 = 1, q2 = 2;
				z1 = Mod.mul_add(z1, z1, C1), z2 = Mod.mul_add(z2, z2, C2);
				for (u64 i = 0; i < M; ++i) {
					u64 t1 = x1 - z1, t2 = x2 - z2;
					z1 = Mod.mul_add(z1, z1, C1), z2 = Mod.mul_add(z2, z2, C2);
					q1 = Mod.mul(q1, t1), q2 = Mod.mul(q2, t2);
				}
				q1 = Mod.mul(q1, x1 - z1), q2 = Mod.mul(q2, x2 - z2);
				const u64 q3 = Mod.mul(q1, q2), g3 = std::gcd(n, q3);
				if (g3 == 1) {
					continue;
				}
				if (g3 != n) {
					ans = g3;
					return;
				}
				const u64 g1 = std::gcd(n, q1);
				const u64 g2 = std::gcd(n, q2);
				const u64 C = g1 != 1 ? C1 : C2;
				const u64 x = g1 != 1 ? x1 : x2;
				u64 z = g1 != 1 ? y1 : y2;
				u64 g = g1 != 1 ? g1 : g2;
				if (g == n) {
					do {
						z = Mod.mul_add(z, z, C);
						g = std::gcd(n, x - z);
					} while (g == 1);
				}
				if (g != n) {
					ans = g;
					return;
				}
				Z1 += 2, Z2 += 2;
				return;
			}
		}
  };
  do {
  	find();
  } while (!ans);
  return ans;
}
template<bool sorted>
vector<pair<u64, u32>> Factorize(u64 n) {
	vector<pair<u64, u32>> ans;
	if (n % 2 == 0) {
		u32 z = __builtin_ctzll(n);
		ans.push_back({2ULL, z}), n >>= z;
	}
	auto upd = [&](const u64& x) {
		for (auto &[p, c] : ans) {
			if (x == p) {
				++ c;
				return;
			}
		}
		ans.push_back({x, 1});
	};
	auto DFS = [&](auto &&self, const u64 &n) -> void {
		if (Is_Prime(n)) {
			return upd(n);
		}
		u64 d = Find_Prime_Factor(n);
		self(self, d), self(self, n / d);
	};
	if (n > 1) {
		DFS(DFS, n); 
	}
	if constexpr (sorted) {
		sort(ans.begin(), ans.end());
	}
	return ans;
}
```