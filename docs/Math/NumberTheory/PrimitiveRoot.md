求原根，需前置 `Factorize()`。

```cpp
u64 primitive_root(u64 P) {
  if (P == 2) return 1;
  if (P == 167772161) return 3;
  if (P == 469762049) return 3;
  if (P == 754974721) return 11;
  if (P == 998244353) return 3;
  auto f = Factorize<false>(P - 1);
  const Montgomery<u64, u128> Mod(P);
  u64 g = Mod.raw(2);
  while (true) {
    if ([&]() {
      for (auto [p, e] : f) {
        if (Mod.same(Mod.pow(g, (P - 1) / p), Mod.one())) {
          return false;
        }
      }
      return true;
    }()) return Mod.val(g);
    g = Mod.inc(g);
  }
}
```