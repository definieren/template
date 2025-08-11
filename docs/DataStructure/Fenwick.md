需要重载运算符 `+=`。

```cpp
template<class T> struct Fenwick {
  private:
    int n, id; vector<T> bit; vector<int> vis;
    
    inline constexpr int hole(int k) {
      return k + (k >> 10);
    }
  
    constexpr void upd(int x, T k) {
      for (x ++; x <= n; x += x & - x) {
//      for (x ++; x; x -= x & - x) {
        if (vis[hole(x) - 1] ^ id) vis[hole(x) - 1] = id, bit[hole(x) - 1] = T{};
        bit[hole(x) - 1] += k;
      }
      return;
    }
    constexpr T qry(int x) {
      T ans{};
      for (x ++; x; x -= x & - x)
//      for (x ++; x <= n; x += x & - x)
        if (vis[hole(x) - 1] == id) ans += bit[hole(x) - 1];
      return ans;
    }
  public:
    constexpr Fenwick(): n(0), id(0), bit{}, vis{} {}
    constexpr Fenwick(int _n): n(_n) {
      id = 0; vis.assign(hole(n), 0);
      bit.assign(hole(n), T{});
    }
    constexpr void Clear() { ++ id; return; }
    
//    Single point addition, prefix / suffix summation.
    constexpr void Update(int x, T k) { return upd(x, k); }
    constexpr T Query(int x) { return qry(x); }
    
//    Single point addition, interval summation.
//    constexpr void Update(int x, T k) { return upd(x, k); }
//    constexpr T Query(int l, int r) { return qry(r) - qry(l - 1); }

//    Interval addition, single point summation.
//    constexpr void Update(int l, int r, T k) { upd(l, k), upd(r + 1, - k); return; }
//    constexpr T Query(int x) { return qry(x); }
};
```