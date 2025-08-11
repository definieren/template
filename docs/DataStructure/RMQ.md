$O(n) - O(1)$ RMQ。

```cpp
template<class T, class Cmp = less<T>>
struct RMQ {
  static constexpr unsigned B = 64;
  const Cmp cmp = Cmp();
  unsigned n;
  vector<vector<T>> st;
  vector<T> pre, suf, A;
  vector<unsigned long long> stk;
  
  RMQ() {}
  RMQ(const vector<T> &a) {
    Init(a);
  }
  
  inline void chkmin(T &a, T b) {
    if (cmp(b, a)) {
      a = b;
    }
    return;
  }
  void Init(const vector<T> &a) {
    n = a.size();
    pre = suf = A = a;
    stk.resize(n);
    if (!n) {
      return;
    }
    const int Bn = (n - 1) / B + 1;
    const int LG = __lg(Bn);
    st.assign(LG + 1, vector<T>(Bn));
    for (int i = 0; i < Bn; i ++) {
      st[0][i] = a[i * B];
      for (unsigned j = 1; j < B && i * B + j < n; j ++) {
        chkmin(st[0][i], a[i * B + j]);
      }
    }
    for (unsigned i = 1; i < n; i ++) {
      if (i % B) {
        chkmin(pre[i], pre[i - 1]);
      }
    }
    for (int i = n - 2; i >= 0; i --) {
      if (i % B != B - 1) {
        chkmin(suf[i], suf[i + 1]);
      }
    }
    for (int i = 0; i < LG; i ++) {
      for (int j = 0; j + (2 << i) <= Bn; j ++) {
        st[i + 1][j] = min(st[i][j], st[i][j + (1 << i)], cmp);
      }
    }
    for (int i = 0; i < Bn; i ++) {
      const int l = i * B;
      const int r = min(n, l + B);
      unsigned long long msk = 0;
      for (int j = l; j < r; j ++) {
        while (msk && cmp(a[j], a[__lg(msk) + l])) {
          msk ^= 1ULL << __lg(msk);
        }
        msk |= 1ULL << (j - l);
        stk[j] = msk;
      }
    }
    return;
  }
  inline int Query(int l, int r) {
    if (l / B != (r - 1) / B) {
      T ans = min(suf[l], pre[r - 1], cmp);
      l = l / B + 1, r = r / B;
      if (l < r) {
        int k = __lg(r - l);
        chkmin(ans, min(st[k][l], st[k][r - (1 << k)], cmp));
      }
      return ans;
    } else {
      int x = B * (l / B);
      return A[__builtin_ctzll(stk[r - 1] >> (l - x)) + l];
    }
  }
};
```