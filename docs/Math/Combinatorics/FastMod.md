用来假装自己没贺板子。

```cpp
constexpr u32 add(const u32& x, const u32& y, const u32& M = Mod) {
  return (x + y >= M) ? (x + y - M) : (x + y);
}
constexpr u32 sub(const u32& x, const u32& y, const u32& M = Mod) {
  return (x - y >= M) ? (x - y + M) : (x - y);
}
constexpr u32 mul(const u32& x, const u32& y, const u32& M = Mod) {
  return static_cast<u32>(static_cast<u64>(x) * y % static_cast<u64>(M));
}
constexpr u32& cadd(u32& x, const u32& y, const u32& M = Mod) {
  return x = add(x, y, M);
}
constexpr u32& csub(u32& x, const u32& y, const u32& M = Mod) {
  return x = sub(x, y, M);
}
constexpr u32& cmul(u32& x, const u32& y, const u32& M = Mod) {
  return x = mul(x, y, M);
}
```