# 圆方树

需前置 `GraphBase.md`。

`block_cut(G)`：返回图 $G$ 的圆方树，方点编号为 $[n, n + \# \text{block})$，圆点编号为 $[0, n)$。

```cpp
template<class Graph>
graph<> block_cut(const Graph& G) {
  assert(!G.is_directed);
  assert(G.is_prepared());
  
  const int n = G.n;
  std::vector<int> low(n), dfn(n, -1), stk;
  graph<> bct(2 * n);
  stk.reserve(n);
  int cur = n, dfc = 0;
  
  auto dfs = [&](auto&& self, int u) -> void {
    low[u] = dfn[u] = dfc ++;
    stk.push_back(u);
    for (auto e : G[u]) {
      const int v = e.v;
      if (dfn[v] == -1) {
        self(self, v);
        low[u] = std::min(low[u], low[v]);
        if (low[v] == dfn[u]) {
          for (int x = -1; x != v; stk.pop_back()) {
            x = stk.back();
            bct.addEdge(cur, x);
          }
          bct.addEdge(cur, u);
          cur ++;
        }
      } else {
        low[u] = std::min(low[u], dfn[v]);
      }
    }
  };
  
  for (int u = 0; u < n; u ++) {
    stk.clear();
    if (dfn[u] == -1) {
      dfs(dfs, u);
      if (dfn[u] + 1 == dfc) {
        bct.addEdge(cur, u);
        cur ++;
      }
    }
  }
  bct.n = cur, bct.build();
  
  return bct;
}
```