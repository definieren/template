```cpp
template<class U0, class U1, class S0, U0 P>
struct Static_Modint {
private:
  static_assert((P >> (sizeof(U0) * 8 - 1)) == 0, "'Mod' must less than max(U0)/2");
  static constexpr U0 Mod = P;
  U0 x;
  template<class T>
  static constexpr unsigned SafeMod(T x) {
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
  constexpr Static_Modint(): x(static_cast<U0>(0)) {}
  template<class T>
  constexpr Static_Modint(T _x): x(SafeMod(_x)) {}
  static constexpr Static_Modint raw(U0 _x) noexcept {
    Static_Modint x;
    return x.x = _x, x;
  }
  static constexpr U0 GetMod() {
    return Mod;
  }
  template<class T>
  explicit constexpr operator T() const {
    return static_cast<T>(x);
  }
  constexpr Static_Modint &operator += (const Static_Modint &rhs) {
    x = ((x += rhs.x) >= Mod) ? (x - Mod) : x;
    return *this;
  }
  constexpr Static_Modint &operator -= (const Static_Modint &rhs) {
    x = ((x -= rhs.x) >= Mod) ? (x + Mod) : x;
    return *this;
  }
  constexpr Static_Modint &operator *= (const Static_Modint &rhs) {
    x = (static_cast<U1>(x) * rhs.x) % Mod;
    return *this;
  }
  constexpr Static_Modint &operator /= (const Static_Modint &rhs) {
    return (*this *= rhs.inv());
  }
  friend constexpr Static_Modint fma(const Static_Modint &a, const Static_Modint &b, const Static_Modint &c) {
    return raw((static_cast<U1>(a.x) * b.x + c.x) % Mod);
  }
  friend constexpr Static_Modint fam(const Static_Modint &a, const Static_Modint &b, const Static_Modint &c) {
    return raw((a.x + static_cast<U1>(b.x) * c.x) % Mod);
  }
  friend constexpr Static_Modint fms(const Static_Modint &a, const Static_Modint &b, const Static_Modint &c) {
    return raw((static_cast<U1>(a.x) * b.x + Mod - c.x) % Mod);
  }
  friend constexpr Static_Modint fsm(const Static_Modint &a, const Static_Modint &b, const Static_Modint &c) {
    return raw((a.x + static_cast<U1>(Mod - b.x) * c.x) % Mod);
  }
  constexpr Static_Modint div_2() const {
    return raw(((x & 1) ? (x + Mod) : x) >> 1);
  }
  constexpr Static_Modint inv() const {
    U0 a = Mod, b = x; S0 y = 0, z = 1;
    while (b) {
      const U0 q = a / b;
      const U0 c = a - q * b;
      a = b, b = c;
      const S0 w = y - static_cast<S0>(q) * z;
      y = z, z = w;
    }
    return raw(y < 0 ? y + Mod : y);
  }
  friend constexpr Static_Modint operator + (const Static_Modint &x) {
    return x;
  }
  friend constexpr Static_Modint operator - (Static_Modint x) {
    x.x = x.x ? (Mod - x.x) : 0U;
    return x;
  }
  constexpr Static_Modint &operator ++ () {
    x = (x + 1 == Mod) ? 0U : (x + 1);
    return *this;
  }
  constexpr Static_Modint &operator -- () {
    x = (x == 0U) ? (Mod - 1) : (x - 1);
    return *this;
  }
  constexpr Static_Modint operator ++ (int) {
    Static_Modint tmp = (*this);
    return ++ (*this), tmp;
  }
  constexpr Static_Modint operator -- (int) {
    Static_Modint tmp = (*this);
    return -- (*this), tmp;
  }
  friend constexpr Static_Modint operator + (Static_Modint x, const Static_Modint &y) {
    return x += y;
  }
  friend constexpr Static_Modint operator - (Static_Modint x, const Static_Modint &y) {
    return x -= y;
  }
  friend constexpr Static_Modint operator * (Static_Modint x, const Static_Modint &y) {
    return x *= y;
  }
  friend constexpr Static_Modint operator / (Static_Modint x, const Static_Modint &y) {
    return x /= y;
  }
  constexpr Static_Modint Pow(long long y) const {
    if (y < 0) return inv().Pow(- y);
    Static_Modint x = *this, ans;
    ans.x = static_cast<U0>(1);
    for (; y; y >>= 1, x *= x)
      if (y & 1) ans *= x;
    return ans;
  }
  friend constexpr std::istream& operator >> (std::istream& is, Static_Modint &x) {
    return is >> x.x;
  }
  friend constexpr std::ostream& operator << (std::ostream& os, const Static_Modint &x) {
    return os << x.x;
  }
  friend constexpr bool operator == (const Static_Modint &x, const Static_Modint &y) {
    return x.x == y.x;
  }
  friend constexpr bool operator != (const Static_Modint &x, const Static_Modint &y) {
    return x.x != y.x;
  }
  friend constexpr bool operator <= (const Static_Modint &x, const Static_Modint &y) {
    return x.x <= y.x;
  }
  friend constexpr bool operator >= (const Static_Modint &x, const Static_Modint &y) {
    return x.x >= y.x;
  }
  friend constexpr bool operator < (const Static_Modint &x, const Static_Modint &y) {
    return x.x < y.x;
  }
  friend constexpr bool operator > (const Static_Modint &x, const Static_Modint &y) {
    return x.x > y.x;
  }
};
template<u32 P>
using sm32 = Static_Modint<u32, u64, int, P>;
template<u64 P>
using sm64 = Static_Modint<u64, u128, i64, P>;
using Z = sm32<998244353U>;
```