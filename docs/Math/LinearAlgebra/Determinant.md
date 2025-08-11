# 行列式

```cpp
mint Det(vector<vector<mint>> a) {
  const int n = static_cast<int>(a.size());
  bool typ = false; mint det = 1;
  for (int i = 0; i < n; i ++) {
    int pos = - 1;
    for (int j = i; j < n; j ++) {
      if (a[j][i]) {
        pos = j; break;
      }
    }
    if (!~pos) {
      return 0;
    }
    if (pos != i) {
      typ ^= 1;
    }
    swap(a[i], a[pos]);
    for (int j = i + 1; j < n; j ++) {
      while (a[j][i]) {
        typ ^= 1;
        int coef = (int)a[i][i] / (int)a[j][i];
        for (int k = i; k < n; k ++) {
          a[i][k] -= a[j][k] * coef, swap(a[i][k], a[j][k]);
        }
      }
    }
    det *= a[i][i];
  }
  return typ ? - det : det;
}
```