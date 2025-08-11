# 自定义 Hash 函数

```cpp
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};
```

# 哈希表

```cpp
template<class Key, class T, class Hash = custom_hash>
struct HashMap {
private:
  u32 rest, msk;
  std::vector<Key> key;
  std::vector<T> val;
  std::vector<bool> used;
  
  u64 H(const Key& x) {
    return Hash{}(x) & msk;
  }
  void rebuild() {
    std::vector<std::pair<Key, T>> d;
    d.reserve(used.size() / 2  - rest);
    for (u32 i = 0; i < used.size(); i ++) {
      if (used[i]) {
        d.push_back({key[i], val[i]});
      }
    }
    reserve(2 * d.size());
    for (auto& [x, y] : d) {
      (*this)[x] = y;
    }
  }
  u32 index(const Key& k) {
    u32 i = H(k);
    while (used[i] && key[i] != k) {
      (++ i) &= msk;
    }
    return i;
  }
public:
  HashMap() {
    reserve(0);
  }
  template<class _InputIterator, class = std::_RequireInputIter<_InputIterator>>
  explicit HashMap(_InputIterator __first, _InputIterator __last) {
    reserve(std::distance(__first, __last));
    for (auto it = __first; it != __last; it ++) {
      (*this)[it -> first] = it -> second;
    }
  }
  explicit HashMap(const std::vector<std::pair<Key, T>> &a) {
    reserve(a.size());
    for (auto i : a) {
      (*this)[i.first] = i.second;
    }
  }
  void reserve(u32 n) {
    u32 n_ = 8;
    while (n_ < n * 2) {
      n_ <<= 1;
    }
    rest = n_ / 2, msk = n_ - 1;
    key.resize(n_);
    val.resize(n_);
    used.assign(n_, 0);
  }
  void clear() {
    used.assign(used.size(), 0);
    rest = (msk + 1) / 2;
  }
  auto size() const {
    return used.size() / 2 - rest;
  }
  bool empty() const {
    return !size();
  }
  bool contains(const Key& k) const {
    u32 i = index(k);
    return used[i] && key[i] == k;
  }
  T& operator [] (const Key& k) {
    if (rest == 0) {
      rebuild();
    }
    const u32 i = index(k);
    if (!used[i]) {
      used[i] = 1;
      key[i] = k;
      val[i] = T{};
      -- rest;
    }
    return val[i];
  }
  T operator [] (const Key& k) const {
    const u32 i = index(k);
    if (!used[i]) {
      return T{};
    }
    return val[i];
  }
  std::vector<std::pair<Key, T>> data() const {
    std::vector<std::pair<Key, T>> d;
    d.reserve(size());
    for (u32 i = 0; i < used.size(); i ++) {
      if (used[i]) {
        d.push_back({key[i], val[i]});
      }
    }
    return d;
  }
  template<class F = void(*)(Key, T&)>
  void for_each(const F& f) {
    for (u32 i = 0; i < used.size(); i ++) {
      if (used[i]) {
        f(key[i], val[i]);
      }
    }
  }
};
```