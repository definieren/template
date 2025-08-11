# SMAWK

完全单调矩阵，求每行最大值位置。

```cpp
// max_i (A)
template<typename T>
vector<int> SMAWK(int n, int m, auto F) {
  vector<int> ans(n, -1);
  auto Solve = [&](auto self, int n, int m, int* const X, int* const Y) -> void {
    if (n <= 2 || m <= 2) {
      for (int i = 0; i < n; ++ i) {
        int x = X[i]; T mx{};
        for (int j = 0; j < m; ++ j) {
          int y = Y[j]; T w = F(x, y);
          if (ans[x] == -1 || w > mx)
            ans[x] = y, mx = w;
        }
      }
      return;
    }
    if (n < m) {
      int k = 0;
      for (int i = 0; i < m; ++ i) {
        int y = Y[i];
        while (k && F(X[k - 1], Y[k - 1]) < F(X[k - 1], y)) -- k;
        if (k < n) Y[k ++] = y;
      }
      m = k;
    }
    auto _X = X + n, _Y = Y + m;
    int _n = 0;
    for (int i = 1; i < n; i += 2)
      _X[_n ++] = X[i];
    for (int i = 0; i < m; ++ i)
      _Y[i] = Y[i];
    self(self, _n, m, _X, _Y);
    int k = 0;
    for (int i = 0; i < n; i += 2) {
      int lim = i + 1 < n ? ans[X[i + 1]] : Y[m - 1];
      T mx{};
      while (true) {
        T w = F(X[i], Y[k]);
        if (ans[X[i]] == -1 || w > mx)
          ans[X[i]] = Y[k], mx = w;
        if (Y[k] == lim) break;
        ++ k;
      }
    }
    return;
  };
  vector<int> X(n << 1), Y(max(m, n << 1));
  for (int i = 0; i < n; ++ i) X[i] = i;
  for (int i = 0; i < m; ++ i) Y[i] = i;
  Solve(Solve, n, m, X.data(), Y.data());
  return ans;
}
```