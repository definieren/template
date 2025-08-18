```cpp
template<class U0, class U1, class S0, U0 P>
struct static_modint {
private:
	static_assert((P >> (sizeof(U0) * 8 - 1)) == 0, "'Mod' must less than max(U0)/2");
	static constexpr U0 Mod = P;
	U0 x;
	template<class T>
	static constexpr U0 safeMod(T x) {
		if constexpr (std::is_unsigned<T>::value) {
			x %= static_cast<T>(Mod);
			return static_cast<U0>(x);
		} else {
			if ((x %= static_cast<T>(Mod)) < 0)
				x += static_cast<T>(Mod);
			return static_cast<U0>(x);
		}
	}
public:
	constexpr static_modint(): x(static_cast<U0>(0)) {}
	template<class T>
	constexpr static_modint(T _x): x(safeMod(_x)) {}
	static constexpr static_modint from_raw(U0 _x) noexcept {
		static_modint x;
		return x.x = _x, x;
	}
	static constexpr U0 getMod() {
		return Mod;
	}
	template<class T>
	explicit constexpr operator T() const {
		return static_cast<T>(x);
	}
	constexpr static_modint &operator += (const static_modint &rhs) {
		x += rhs.x, (x - Mod) >> (sizeof(U0) * 8 - 1) || (x -= Mod);
		return *this;
	}
	constexpr static_modint &operator -= (const static_modint &rhs) {
		x -= rhs.x, x >> (sizeof(U0) * 8 - 1) && (x += Mod);
		return *this;
	}
	constexpr static_modint &operator *= (const static_modint &rhs) {
		x = static_cast<U0>(static_cast<U1>(x) * static_cast<U1>(rhs.x) % static_cast<U1>(Mod));
		return *this;
	}
	constexpr static_modint &operator /= (const static_modint &rhs) {
		return (*this *= rhs.inv());
	}
	friend constexpr static_modint fma(const static_modint &a, const static_modint &b, const static_modint &c) {
		return from_raw((static_cast<U1>(a.x) * b.x + c.x) % Mod);
	}
	friend constexpr static_modint fam(const static_modint &a, const static_modint &b, const static_modint &c) {
		return from_raw((a.x + static_cast<U1>(b.x) * c.x) % Mod);
	}
	friend constexpr static_modint fms(const static_modint &a, const static_modint &b, const static_modint &c) {
		return from_raw((static_cast<U1>(a.x) * b.x + Mod - c.x) % Mod);
	}
	friend constexpr static_modint fsm(const static_modint &a, const static_modint &b, const static_modint &c) {
		return from_raw((a.x + static_cast<U1>(Mod - b.x) * c.x) % Mod);
	}
	constexpr static_modint inv() const {
		U0 a = Mod, b = x; S0 y = 0, z = 1;
		while (b) {
			const U0 q = a / b;
			const U0 c = a - q * b;
			a = b, b = c;
			const S0 w = y - static_cast<S0>(q) * z;
			y = z, z = w;
		}
		return from_raw(y < 0 ? y + Mod : y);
	}
	friend constexpr static_modint inv(const static_modint &x) {
		return x.inv();
	}
	friend constexpr static_modint operator + (const static_modint &x) {
		return x;
	}
	friend constexpr static_modint operator - (static_modint x) {
		x.x = x.x ? (Mod - x.x) : 0U;
		return x;
	}
	constexpr static_modint &operator ++ () {
		return *this += 1;
	}
	constexpr static_modint &operator -- () {
		return *this -= 1;
	}
	constexpr static_modint operator ++ (int) {
		static_modint v = *this;
		return *this += 1, v;
	}
	constexpr static_modint operator -- (int) {
		static_modint v = *this;
		return *this -= 1, v;
	}
	friend constexpr static_modint operator + (static_modint x, const static_modint &y) {
		return x += y;
	}
	friend constexpr static_modint operator - (static_modint x, const static_modint &y) {
		return x -= y;
	}
	friend constexpr static_modint operator * (static_modint x, const static_modint &y) {
		return x *= y;
	}
	friend constexpr static_modint operator / (static_modint x, const static_modint &y) {
		return x /= y;
	}
	template<class T>
	constexpr static_modint pow(T y) const {
		if (y < 0) return inv().pow(- y);
		static_modint x = *this, ans = from_raw(1U);
		for (; y; y >>= 1, x *= x) {
			if (y & 1) {
				ans *= x;
			}
		}
		return ans;
	}
	template<class T>
	friend constexpr static_modint pow(const static_modint &x, T y) {
		return x.pow(y);
	}
	std::pair<bool, static_modint> sqrt() const {
		if (x == 0U) {
			return {true, from_raw(0)};
		}
		if (Mod == 2U) {
			return {true, from_raw(1)};
		}
		if (pow((Mod - 1) / 2) != from_raw(1)) {
			return {false, from_raw(0)};
		}
		static std::mt19937_64 rnd(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_int_distribution<U0> uid(1U, Mod - 1);
		static_modint x, y;
		do {
			x = from_raw(uid(rnd));
			y = x * x - *this;
		} while (y.pow((Mod - 1) / 2) == from_raw(1));
		auto mul = [](std::pair<static_modint, static_modint> &f,
			const std::pair<static_modint, static_modint> &g, const static_modint &h) {
			f = {f.first * g.first + f.second * g.second * h,
					 f.first * g.second + f.second * g.first};
		};
		std::pair<static_modint, static_modint> f{x, 1}, g{1, 0};
		auto exp = (Mod + 1) / 2;
		for (; exp; exp >>= 1, mul(f, f, y)) {
			if (exp & 1) {
				mul(g, f, y);
			}
		}
		return {true, from_raw(std::min(g.first.x, Mod - g.first.x))};
	}
	friend std::pair<bool, static_modint> sqrt(const static_modint &x) {
		return x.sqrt();
	}
	friend constexpr std::istream& operator >> (std::istream& is, static_modint &x) {
		S0 y;
		is >> y, x = y;
		return is;
	}
	friend constexpr std::ostream& operator << (std::ostream& os, const static_modint &x) {
		return os << x.x;
	}
	friend constexpr bool operator == (const static_modint &x, const static_modint &y) {
		return x.x == y.x;
	}
	friend constexpr bool operator != (const static_modint &x, const static_modint &y) {
		return x.x != y.x;
	}
	friend constexpr bool operator <= (const static_modint &x, const static_modint &y) {
		return x.x <= y.x;
	}
	friend constexpr bool operator >= (const static_modint &x, const static_modint &y) {
		return x.x >= y.x;
	}
	friend constexpr bool operator < (const static_modint &x, const static_modint &y) {
		return x.x < y.x;
	}
	friend constexpr bool operator > (const static_modint &x, const static_modint &y) {
		return x.x > y.x;
	}
};
template<u32 P>
using sm32 = static_modint<u32, u64, int, P>;
template<u64 P>
using sm64 = static_modint<u64, u128, i64, P>;
using Z = sm32<998244353U>;
```