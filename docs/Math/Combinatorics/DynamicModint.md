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
	template<class T>
	constexpr U0 trans(T x) const noexcept {
		if (__builtin_expect(static_cast<T>(0) <= x && x < static_cast<T>(n), 1)) {
			return raw(static_cast<U0>(x));
		}
		if constexpr (is_unsigned<T>::value) {
			x %= static_cast<T>(n);
		} else {
			if ((x %= static_cast<T>(n)) < 0) {
				(x += static_cast<T>(n)) %= static_cast<T>(n);
			}
		}
		return Reduce(static_cast<U0>(x), rs);
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

template<class U0, class U1, class S0, unsigned P>
struct Dynamic_Modint {
private:
	static inline Montgomery<U0, U1> Mod = P;
	U0 x;
public:
	constexpr Dynamic_Modint(): x(Mod.zero()) {}
	template<class T>
	constexpr Dynamic_Modint(T _x): x(Mod.trans(_x)) {}
	static constexpr Dynamic_Modint raw(U0 _x) {
		Dynamic_Modint x;
		return x.x = Mod.raw(_x), x;
	}
	static constexpr U0 GetMod() {
		return Mod.GetMod();
	}
	static constexpr void SetMod(U0 Mod_) {
		return Mod.SetMod(Mod_);
	}
	template<class T>
	explicit constexpr operator T() const {
		return static_cast<T>(Mod.val(x));
	}
	constexpr Dynamic_Modint &operator += (const Dynamic_Modint &rhs) {
		x = Mod.add(x, rhs.x);
		return *this;
	}
	constexpr Dynamic_Modint &operator -= (const Dynamic_Modint &rhs) {
		x = Mod.sub(x, rhs.x);
		return *this;
	}
	constexpr Dynamic_Modint &operator *= (const Dynamic_Modint &rhs) {
		x = Mod.mul(x, rhs.x);
		return *this;
	}
	constexpr Dynamic_Modint &operator /= (const Dynamic_Modint &rhs) {
		return (*this *= rhs.inv());
	}
	constexpr Dynamic_Modint inv() const {
		U0 a = GetMod(), b = Mod.val(x);
		S0 y = 0, z = 1;
		while (b) {
			const U0 q = a / b;
			const U0 c = a - q * b;
			a = b, b = c;
			const S0 w = y - static_cast<S0>(q) * z;
			y = z, z = w;
		}
		return raw((y < 0) ? (y + GetMod()) : y);
	}
	friend constexpr Dynamic_Modint operator + (const Dynamic_Modint &x) {
		return x;
	}
	friend constexpr Dynamic_Modint operator - (Dynamic_Modint x) {
		return x.x = Mod.neg(x.x), x;
	}
	constexpr Dynamic_Modint &operator ++ () {
		x = Mod.inc(x);
		return *this;
	}
	constexpr Dynamic_Modint &operator -- () {
		x = Mod.dec(x);
		return *this;
	}
	constexpr Dynamic_Modint operator ++ (int) {
		Dynamic_Modint t = (*this);
		return ++ (*this), t;
	}
	constexpr Dynamic_Modint operator -- (int) {
		Dynamic_Modint t = (*this);
		return -- (*this), t;
	}
	friend constexpr Dynamic_Modint operator + (Dynamic_Modint x, const Dynamic_Modint &y) {
		return x += y;
	}
	friend constexpr Dynamic_Modint operator - (Dynamic_Modint x, const Dynamic_Modint &y) {
		return x -= y;
	}
	friend constexpr Dynamic_Modint operator * (Dynamic_Modint x, const Dynamic_Modint &y) {
		return x *= y;
	}
	friend constexpr Dynamic_Modint operator / (Dynamic_Modint x, const Dynamic_Modint &y) {
		return x /= y;
	}
	constexpr Dynamic_Modint Pow(long long y) const {
		if (y < 0) return inv().Pow(- y);
		Dynamic_Modint x = *this, ans;
		ans.x = Mod.one();
		for (; y; y >>= 1, x *= x)
			if (y & 1) ans *= x;
		return ans;
	}
	friend ostream& operator << (ostream& os, const Dynamic_Modint &x) {
		return os << Mod.val(x.x);
	}
	friend constexpr bool operator == (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return Mod.same(x.x, y.x);
	}
	friend constexpr bool operator != (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return !(x == y);
	}
	friend constexpr bool operator <= (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return Mod.val(x.x) <= Mod.val(y.x);
	}
	friend constexpr bool operator >= (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return Mod.val(x.x) >= Mod.val(y.x);
	}
	friend constexpr bool operator < (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return Mod.val(x.x) < Mod.val(y.x);
	}
	friend constexpr bool operator > (const Dynamic_Modint &x, const Dynamic_Modint &y) {
		return Mod.val(x.x) > Mod.val(y.x);
	}
};
template<u32 P>
using dm32 = Dynamic_Modint<u32, u64, int, P>;
template<u64 P>
using dm64 = Dynamic_Modint<u64, u128, i64, P>;
using Z = dm32<MOD>;
```