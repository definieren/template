# 连续点值

## 拉格朗日插值求多项式系数

```cpp
// x = 0, 1, 2, ..., n - 1
vector<Z> Lagrange(vector<Z> y) {
  if (y.empty()) return y;
  const int n = (int)y.size();
  vector<Z> f(n + 1), a(n), coef, g(n);
  
  for (int i = 0; i < n; i ++) {
    Z res = y[i] * comb.ifac(i) * comb.ifac(n - i - 1);
    ((n - i - 1) & 1) ? a[i] = - res : a[i] = res;
  }
  
  f[0] = 1;
  for (int i = 0; i < n; i ++) {
    f[i + 1] = 0;
    for (int j = i + 1; j; j --)
      f[j] = f[j - 1] - i * f[j];
    f[0] = - f[0] * i;
  }
  
  coef = f, coef.erase(coef.begin());
  for (auto &x : coef) x *= a[0];
  for (int i = 1; i < n; i ++) {
    g[0] = - comb.inv(i) * f[0];
    for (int j = 1; j < n; j ++)
      g[j] = comb.inv(i) * (g[j - 1] - f[j]);
    for (int j = 0; j < n; j ++)
      coef[j] += a[i] * g[j];
  }
  
  return coef;
}
```
