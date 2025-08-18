```cpp
template<class U0, class U1, class S0, U0 P_>
struct montgomery_modint {
private:
	static_assert((P_ & 1), "Mod must be odd");
	static_assert((P_ >> (sizeof(U0) * 8 - 2)) == 0, "'Mod' must less than max(U0)/4");
	template<class T>
	static constexpr U0 safeMod(T x) {
		if constexpr (std::is_unsigned<T>::value) {
			x %= static_cast<T>(P);
			return static_cast<U0>(x);
		} else {
			if ((x %= static_cast<T>(P)) < 0)
				x += static_cast<T>(P);
			return static_cast<U0>(x);
		}
	}
	static constexpr U0 uinv(U0 x) {
		U0 y = x;
		for (int i = sizeof(U0); i; i --) {
			y *= 2 - x * y;
		}
		return y;
	}
	static constexpr U0 reduce(U1 x) {
		return static_cast<U0>((x + static_cast<U0>(x) * R * static_cast<U1>(P)) >> (sizeof(U0) * 8));
	}
	static constexpr U0 P = P_, P2 = P << 1, R = -uinv(P), R2 = (-static_cast<U1>(P)) % P;
	U0 x;
public:
	constexpr montgomery_modint(): x(static_cast<U0>(0)) {}
	template<class T>
	constexpr montgomery_modint(T _x): x(reduce(static_cast<U1>(R2) * safeMod(_x))) {}
	template<int o = 0>
	static constexpr montgomery_modint from_raw(U0 _x) noexcept {
		if constexpr (!o) {
			montgomery_modint x;
			return x.x = reduce(static_cast<U1>(R2) * _x), x;
		} else {
			montgomery_modint x;
			return x.x = _x, x;
		}
	}
	static constexpr U0 getMod() {
		return P;
	}
	constexpr U0 val() const {
		U0 v = reduce(x);
		return (v - P) >> (sizeof(U0) * 8 - 1) ? v : (v - P);
	}
	template<class T>
	explicit constexpr operator T() const {
		return static_cast<T>(val());
	}
	constexpr montgomery_modint &operator += (const montgomery_modint &rhs) {
		x += rhs.x, (x - P2) >> (sizeof(U0) * 8 - 1) || (x -= P2);
		return *this;
	}
	constexpr montgomery_modint &operator -= (const montgomery_modint &rhs) {
		x -= rhs.x, x >> (sizeof(U0) * 8 - 1) && (x += P2);
		return *this;
	}
	constexpr montgomery_modint &operator *= (const montgomery_modint &rhs) {
		x = reduce(static_cast<U1>(x) * rhs.x);
		return *this;
	}
	constexpr montgomery_modint &operator /= (const montgomery_modint &rhs) {
		return (*this *= rhs.inv());
	}
	friend constexpr montgomery_modint fma(const montgomery_modint &a, const montgomery_modint &b, const montgomery_modint &c) {
		return from_raw<1>(reduce(static_cast<U1>(a.x) * b.x + c.x));
	}
	friend constexpr montgomery_modint fam(const montgomery_modint &a, const montgomery_modint &b, const montgomery_modint &c) {
		return from_raw<1>(reduce(a.x + static_cast<U1>(b.x) * c.x));
	}
	friend constexpr montgomery_modint fms(const montgomery_modint &a, const montgomery_modint &b, const montgomery_modint &c) {
		return from_raw<1>(reduce(static_cast<U1>(a.x) * b.x + P2 - c.x));
	}
	friend constexpr montgomery_modint fsm(const montgomery_modint &a, const montgomery_modint &b, const montgomery_modint &c) {
		return from_raw<1>(reduce(a.x + static_cast<U1>(P2 - b.x) * c.x));
	}
	constexpr montgomery_modint inv() const {
		return pow(P - 2);
	}
	friend constexpr montgomery_modint inv(const montgomery_modint &x) {
		return x.inv();
	}
	friend constexpr montgomery_modint operator + (const montgomery_modint &x) {
		return x;
	}
	friend constexpr montgomery_modint operator - (montgomery_modint x) {
		x.x = x.x ? (P2 - x.x) : 0U;
		return x;
	}
	constexpr montgomery_modint &operator ++ () {
		return *this += 1;
	}
	constexpr montgomery_modint &operator -- () {
		return *this -= 1;
	}
	constexpr montgomery_modint operator ++ (int) {
		montgomery_modint v = *this;
		return *this += 1, v;
	}
	constexpr montgomery_modint operator -- (int) {
		montgomery_modint v = *this;
		return *this -= 1, v;
	}
	friend constexpr montgomery_modint operator + (montgomery_modint x, const montgomery_modint &y) {
		return x += y;
	}
	friend constexpr montgomery_modint operator - (montgomery_modint x, const montgomery_modint &y) {
		return x -= y;
	}
	friend constexpr montgomery_modint operator * (montgomery_modint x, const montgomery_modint &y) {
		return x *= y;
	}
	friend constexpr montgomery_modint operator / (montgomery_modint x, const montgomery_modint &y) {
		return x /= y;
	}
	template<class T>
	constexpr montgomery_modint pow(T y) const {
		if (y < 0) return inv().pow(- y);
		montgomery_modint x = *this, ans = from_raw(1U);
		for (; y; y >>= 1, x *= x) {
			if (y & 1) {
				ans *= x;
			}
		}
		return ans;
	}
	template<class T>
	friend constexpr montgomery_modint pow(const montgomery_modint &x, T y) {
		return x.pow(y);
	}
	std::pair<bool, montgomery_modint> sqrt() const {
		static constexpr montgomery_modint one = from_raw(1);
		static constexpr montgomery_modint zero = from_raw(0);
		if (!val()) {
			return {true, zero};
		}
		if (pow((P - 1) / 2) != one) {
			return {false, zero};
		}
		static std::mt19937_64 rnd(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_int_distribution<U0> uid(1U, P - 1);
		montgomery_modint x, y;
		do {
			x = from_raw(uid(rnd));
			y = x * x - *this;
		} while (y.pow((P - 1) / 2) == one);
		auto mul = [](std::pair<montgomery_modint, montgomery_modint> &f,
			const std::pair<montgomery_modint, montgomery_modint> &g, const montgomery_modint &h) {
			f = {f.first * g.first + f.second * g.second * h,
					 f.first * g.second + f.second * g.first};
		};
		std::pair<montgomery_modint, montgomery_modint> f{x, one}, g{one, zero};
		auto exp = (P + 1) / 2;
		for (; exp; exp >>= 1, mul(f, f, y)) {
			if (exp & 1) {
				mul(g, f, y);
			}
		}
		return {true, std::min(g.first, -g.first)};
	}
	friend std::pair<bool, montgomery_modint> sqrt(const montgomery_modint &x) {
		return x.sqrt();
	}
	friend constexpr std::istream& operator >> (std::istream& is, montgomery_modint &x) {
		S0 y;
		is >> y, x = y;
		return is;
	}
	friend constexpr std::ostream& operator << (std::ostream& os, const montgomery_modint &x) {
		return os << x.val();
	}
	friend constexpr bool operator == (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() == y.val();
	}
	friend constexpr bool operator != (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() != y.val();
	}
	friend constexpr bool operator <= (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() <= y.val();
	}
	friend constexpr bool operator >= (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() >= y.val();
	}
	friend constexpr bool operator < (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() < y.val();
	}
	friend constexpr bool operator > (const montgomery_modint &x, const montgomery_modint &y) {
		return x.val() > y.val();
	}
};
template<u32 P>
using mm32 = montgomery_modint<u32, u64, int, P>;
template<u64 P>
using mm64 = montgomery_modint<u64, u128, i64, P>;
using Z = mm32<998244353U>;
```