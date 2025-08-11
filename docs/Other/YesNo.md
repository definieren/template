`cout` 版本：

```cpp
bool YES(bool f = true) {
  std::cout << (f ? "YES" : "NO") << '\n';
  return f;
}
bool NO(bool f = true) {
  return !YES(!f);
}
bool Yes(bool f = true) {
  std::cout << (f ? "Yes" : "No") << '\n';
  return f;
}
bool No(bool f = true) {
  return !Yes(!f);
}
bool yes(bool f = true) {
  std::cout << (f ? "yes" : "no") << '\n';
  return f;
}
bool no(bool f = true) {
  return !yes(!f);
}
```


快写版本：

```cpp
void YES(bool f = true) {
  f ? Puts("YES") : Puts("NO");
  return;
}
void NO(bool f = true) {
  YES(!f);
  return;
}
void Yes(bool f = true) {
  f ? Puts("Yes") : Puts("No");
  return;
}
void No(bool f = true) {
  Yes(!f);
  return;
}
void yes(bool f = true) {
  f ? Puts("yes") : Puts("no");
  return;
}
void no(bool f = true) {
  yes(!f);
  return;
}
```