# 欧拉路径

需前置 `GraphBase`、`HashMap`。

`eulerWalk(G, s)`：求图 $G$ 以 $s$ 为起点的欧拉路径（如果未给定 $s$ 则为自选起点）。

返回值为二元组 $(vs, es)$ 分别表示欧拉路径经过的点、边，返回空则表示不存在欧拉路径。

时间复杂度 $O(n + m)$。


```cpp
template<class Graph, bool allow_use_twice = false>
std::vector<int> vs_to_es(const Graph& G, const std::vector<int>& vs) {
  if (vs.empty()) {
    return {};
  }
  
  const int n = G.n, m = G.m;
  HashMap<u64, int> mp;
  std::vector<int> nxt(m);
  mp.reserve(m);
  
  auto get = [&](int u, int v) -> u64 {
    if (!G.is_directed && u > v) {
      std::swap(u, v);
    }
    return static_cast<u64>(u) * n + v;
  };
  
  for (int i = 0; i < m; i ++) {
    const u64 id = get(G.edges[i].u, G.edges[i].v);
    int i_ = -1;
    if (mp.contains(id)) {
      i_ = mp[id];
    }
    nxt[i] = i_, mp[id] = i;
  }
  
  const int V = vs.size();
  std::vector<int> es(V - 1);
  for (int i = 0; i + 1 < V; i ++) {
    const u64 id = get(vs[i], vs[i + 1]);
    int eid = -1;
    if (mp.contains(id)) {
      eid = mp[id];
    }
    assert(eid != -1);
    es[i] = eid;
    if constexpr (!allow_use_twice) {
      mp[id] = nxt[eid];
    }
  }
  
  return es;
}

template<class Graph>
std::pair<std::vector<int>, std::vector<int>> eulerWalk(const Graph& G, int s = -1) {
  const int n = G.n, m = G.m;
  assert(G.is_prepared());
  assert(n);
  
  if (s == -1) {
    std::vector<int> deg(n);
    for (auto&& e : G.edges) {
      if constexpr (G.is_directed) {
        ++ deg[e.u], -- deg[e.v];
      } else {
        ++ deg[e.u], ++ deg[e.v];
      }
    }
    if constexpr (G.is_directed) {
      s = std::max_element(deg.begin(), deg.end()) - deg.begin();
      if (deg[s] == 0) {
        s = (!m ? 0 : G.edges.front().u);
      }
    } else {
      s = [&]() {
        for (int u = 0; u < n; u ++) {
          if (deg[u] & 1) {
            return u;
          }
        }
        return (!m ? 0 : G.edges.front().u);
      }();
    }
  }
  
  if (m == 0) {
    return {{s}, {}};
  }
  
  std::vector<int> D(n), head = G.head, stk{s}, vs;
  std::vector<bool> vis(m);
  stk.reserve(m + 1), vs.reserve(m + 1);
  ++ D[s];
  while (stk.size()) {
    const int u = stk.back();
    if (head[u] == G.head[u + 1]) {
      vs.push_back(u);
      stk.pop_back();
      continue;
    }
    const auto e = G.csr_edges[head[u] ++];
    const int v = e.v, id = e.id;
    if (!vis[id]) {
      -- D[u], ++ D[v];
      vis[id] = true;
      stk.push_back(v);
    }
  }
  
  for (auto x : D) {
    if (x < 0) {
      return {{}, {}};
    }
  }
  if (static_cast<int>(vs.size()) != m + 1) {
    return {{}, {}};
  }
  std::reverse(vs.begin(), vs.end());
  return {vs, vs_to_es(G, vs)};
}
```