基本用法和 `std::bitset` 一样，区别是修改用 `.set()`，查询用 `.get()`。

```cpp
template<class T>
struct BitSet {
  static_assert(std::is_unsigned_v<T> || std::is_same_v<T, __uint128_t>,
    "The type of bitset must be unsigned");
  
  private:
    static constexpr size_t B = sizeof(T) * 8, lg = std::__lg(B);
    static_assert((B & -B) == B);
    std::vector<T> num;
    
    static constexpr size_t ctz(T x) {
      if constexpr (B <= 32) return __builtin_ctz(x);
      else if constexpr (B == 64) return __builtin_ctzll(x);
      else {
        const unsigned long long lo = x & (~0ULL), hi = x >> 64;
        if (lo) return __builtin_ctzll(lo);
        else return 64 + __builtin_ctzll(hi);
      }
    }
    static constexpr size_t ppc(T x) {
      if constexpr (B <= 32) return __builtin_popcount(x);
      else if constexpr (B == 64) return __builtin_popcountll(x);
      else return __builtin_popcountll(x & (~0ULL))
        + __builtin_popcountll(x >> 64);
    }
    static constexpr size_t which_block(size_t x) {
      if constexpr (B == 8) return x >> 3;
      else if constexpr (B == 16) return x >> 4;
      else if constexpr (B == 32) return x >> 5;
      else if constexpr (B == 64) return x >> 6;
      else return x >> 7;
    }
    static constexpr size_t which_bit(size_t x) {
      if constexpr (B == 8) return x & 7;
      else if constexpr (B == 16) return x & 15;
      else if constexpr (B == 32) return x & 31;
      else if constexpr (B == 64) return x & 63;
      else return x & 127;
    }
  
  public:
    constexpr BitSet(): num{} {}
    explicit constexpr BitSet(const std::vector<T> &vec): num(vec) {}
    explicit constexpr BitSet(const size_t &sz): num(which_block(sz + B - 1)) {}
    explicit constexpr BitSet(const std::initializer_list<T> &a): num(a) {}
    template<class _InputIterator, class = std::_RequireInputIter<_InputIterator>>
    explicit constexpr BitSet(_InputIterator __first, _InputIterator __last): num(__first, __last) {}
    template<class F = bool(*)(int)>
    explicit constexpr BitSet(const size_t &n, F f) {
      resize(n);
      for (size_t i = 0; i < n; i ++) set(i, f(i));
    }
    
    constexpr size_t size() const {
      return num.size() * B;
    }
    constexpr void resize(const size_t &sz) {
      num.resize(which_block(sz + B - 1));
      return;
    }
    
    constexpr BitSet& operator <<= (const size_t &len) {
      const size_t n = num.size();
      if (__builtin_expect(len < size(), 1)) {
        if (__builtin_expect(len != 0, 1)) {
          const size_t a = which_block(len), b = which_bit(len);
          if (b == 0) {
            for (size_t i = n - 1; i >= a; i --)
              num[i] = num[i - a];
          } else {
            const size_t _b = B - b;
            for (size_t i = n - 1; i > a; i --)
              num[i] = (num[i - a] << b) | (num[i - a - 1] >> _b);
            num[a] = num[0] << b;
          }
          std::fill(num.begin(), num.begin() + a, static_cast<T>(0));
        }
      } else reset();
      return *this;
    }
    constexpr BitSet& operator >>= (const size_t &len) {
      const size_t n = num.size();
      if (__builtin_expect(len < size(), 1)) {
        if (__builtin_expect(len != 0, 1)) {
          const size_t a = which_block(len), b = which_bit(len), c = n - a - 1;
          if (b == 0) {
            for (size_t i = 0; i <= c; i ++)
              num[i] = num[i + a];
          } else {
            const size_t _b = B - b;
            for (size_t i = 0; i < c; i ++)
              num[i] = (num[i + a] >> b) | (num[i + a + 1] << _b);
            num[c] = num[n - 1] >> b;
          }
          std::fill(num.begin() + c + 1, num.end(), static_cast<T>(0));
        }
      } else reset();
      return *this;
    }
    constexpr BitSet& operator |= (const BitSet &rhs) {
      const size_t n = size() / B, m = rhs.size() / B;
      num.resize(std::max(n, m)); 
      for (size_t i = 0; i < m; i ++) num[i] |= rhs.num[i];
      return *this;
    }
    constexpr BitSet& operator &= (const BitSet &rhs) {
      const size_t n = size() / B, m = rhs.size() / B;
      num.resize(std::max(n, m)); 
      for (size_t i = 0; i < m; i ++) num[i] &= rhs.num[i];
      return *this;
    }
    constexpr BitSet& operator ^= (const BitSet &rhs) {
      const size_t n = size() / B, m = rhs.size() / B;
      num.resize(std::max(n, m)); 
      for (size_t i = 0; i < m; i ++) num[i] ^= rhs.num[i];
      return *this;
    }
    
    friend constexpr BitSet operator << (BitSet lhs, const size_t &rhs) {
      return lhs <<= rhs;
    }
    friend constexpr BitSet operator >> (BitSet lhs, const size_t &rhs) {
      return lhs >>= rhs;
    }
    friend constexpr BitSet operator | (BitSet lhs, const BitSet &rhs) {
      return lhs |= rhs;
    }
    friend constexpr BitSet operator & (BitSet lhs, const BitSet &rhs) {
      return lhs &= rhs;
    }
    friend constexpr BitSet operator ^ (BitSet lhs, const BitSet &rhs) {
      return lhs ^= rhs;
    }
    friend constexpr BitSet operator ~ (BitSet rhs) {
      const size_t n = rhs.size() / B;
      for (size_t i = 0; i < n; i ++) rhs.num[i] = ~rhs.num[i];
      return rhs;
    }
    
    constexpr T operator [] (const size_t &x) const {
      assert(x < num.size());
      return num[x];
    }
    constexpr T& operator [] (const size_t &x) {
      assert(x < num.size());
      return num[x];
    }
    
    friend constexpr bool operator == (const BitSet &lhs, const BitSet &rhs) {
      const size_t n = lhs.size() / B, m = rhs.size() / B;
      if (n != m) return false;
      for (size_t i = 0; i < n; i ++)
        if (lhs[i] != rhs[i]) return false;
      return true;
    }
    friend constexpr bool operator != (const BitSet &lhs, const BitSet &rhs) {
      return !(lhs == rhs);
    }
    
    friend constexpr std::ostream& operator << (std::ostream& os, const BitSet &x) {
      const int n = x.size();
      for (int i = 0; i < n; i ++) os << x.get(i);
      return os;
    }
    
    constexpr unsigned get(const size_t &x) const {
      assert(x < size());
      return (num[which_block(x)] >> (which_bit(x))) & 1;
    }
    constexpr void flip(const size_t &x) {
      assert(x < size());
      num[which_block(x)] ^= (static_cast<T>(1) << (which_bit(x)));
      return;
    }
    constexpr void set(const size_t &x) {
      assert(x < size());
      num[which_block(x)] |= (static_cast<T>(1) << which_bit(x));
      return;
    }
    constexpr void set(const size_t &x, const unsigned &o) {
      assert(x < size());
      o ? set(x) : reset(x);
      return;
    }
    constexpr void set() {
      __builtin_memset(num.data(), 0xff, sizeof(T) * num.size());
      return;
    }
    constexpr void reset(const size_t &x) {
      assert(x < size());
      num[which_block(x)] &= ~(static_cast<T>(1) << which_bit(x));
      return;
    }
    constexpr void reset() {
      __builtin_memset(num.data(), 0, sizeof(T) * num.size());
      return;
    }
    constexpr size_t count() const {
      size_t ret = 0;
      for (auto i : num) ret += ppc(i);
      return ret;
    }
    constexpr bool none() const {
      for (auto i : num)
        if (i != static_cast<T>(0)) return false;
      return true;
    }
    constexpr bool all() const {
      for (auto i : num)
        if (i != (~static_cast<T>(0))) return false;
      return true;
    }
    constexpr bool any() const {
      return !none();
    }
    
    constexpr size_t Find_First() const {
      const size_t n = num.size();
      for (size_t i = 0; i < n; i ++)
        if (num[i] != static_cast<T>(0))
          return i * B + ctz(num[i]);
      return size();
    }
    constexpr size_t Find_Next(size_t x) const {
      const size_t n = num.size();
      if (++ x >= size()) return size();
      size_t i = which_block(x), j = which_bit(x);
      T numi = num[i];
      numi &= (~static_cast<T>(0)) << j;
      if (numi != static_cast<T>(0))
        return i * B + ctz(numi);
      for (++ i; i < n; i ++) {
        T numi = num[i];
        if (numi) return i * B + ctz(num[i]);
      }
      return size();
    }
}; using Bitset = BitSet<unsigned long long>;
```