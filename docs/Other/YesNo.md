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