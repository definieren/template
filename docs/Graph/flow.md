# 最大流

```cpp
template<class Cap>
struct MaxFlow_Graph {
  private:
    struct _Edge {
      int to, rev;
      Cap cap;
    };
    
    int n;
    vector<pair<int, int>> pos;
    vector<vector<_Edge>> G;
  public:
    struct Edge {
      int u, v;
      Cap cap, flow;
    };
    
    MaxFlow_Graph(): n(0) { G.clear(), pos.clear(); return; }
    explicit MaxFlow_Graph(int _n): n(_n)
      { G.assign(n, vector<_Edge>()), pos.clear(); return; }
    
    int Add_Edge(int u, int v, Cap cap) {
      int m = (int)pos.size();
      pos.emplace_back(u, (int)G[u].size());
      int uid = (int)G[u].size(), vid = (int)G[v].size();
      if (u == v) vid ++;
      G[u].emplace_back(v, vid, cap);
      G[v].emplace_back(u, uid, 0);
      return m;
    }
    Edge Get_Edge(int i) {
      auto e = G[pos[i].fir][pos[i].sec], re = G[e.to][e.rev];
      return Edge{pos[i].fir, e.to, e.cap + re.cap, re.cap};
    }
    vector<Edge> Get_All_Edges() {
      int m = (int)pos.size(); vector<Edge> res;
      for (int i = 0; i < m; i ++) res.emplace_back(Get_Edge(i));
      return res;
    }
    void Change_Edge(int i, Cap _cap, Cap _flow) {
      auto &e = G[pos[i].fir][pos[i].sec], &re = G[e.to][e.rev];
      e.cap = _cap - _flow, re.cap = _flow; return;
    }
    Cap Max_Flow(int S, int T) {
      return Max_Flow(S, T, numeric_limits<Cap>::max());
    }
    Cap Max_Flow(int S, int T, Cap lim) {
      Cap flow = 0;
      vector<int> dep(n), cur(n); queue<int> Q;
      auto BFS = [&]() {
        fill(dep.begin(), dep.end(), -1);
        while (Q.size()) Q.pop();
        Q.emplace(S), dep[S] = 0;
        while (Q.size()) {
          int u = Q.front(); Q.pop();
          for (auto e : G[u]) {
            if (e.cap == 0 || ~dep[e.to]) continue;
            dep[e.to] = dep[u] + 1;
            if (e.to == T) return true;
            Q.emplace(e.to);
          }
        }
        return false;
      };
      auto DFS = [&](auto self, int u, Cap flow) {
        if (u == T) return flow;
        Cap used = 0;
        for (int &i = cur[u]; i < (int)G[u].size(); i ++) {
          auto &e = G[u][i], &re = G[e.to][e.rev];
          if (dep[e.to] == dep[u] + 1 && e.cap) {
            Cap d = self(self, e.to, min(flow - used, e.cap));
            if (d > 0) e.cap -= d, re.cap += d, used += d;
            else dep[e.to] = -1;
            if (flow == used) break;
          }
        }
        return used;
      };
      while (flow < lim && BFS()) {
        for (int i = 0; i < n; i ++) cur[i] = 0;
        Cap f = DFS(DFS, S, lim - flow);
        if (!f) break;
        flow += f;
      }
      return flow;
    }
    vector<bool> Min_Cut(int S) {
      vector<bool> vis(n); queue<int> Q;
      Q.emplace(S), vis[S] = true;
      while (Q.size()) {
        int u = Q.front(); Q.pop();
        for (auto e : G[u])
          if (e.cap && !vis[e.to])
            vis[e.to] = true, Q.emplace(e.to);
      }
      return vis;
    }
};
```

# 最小费用最大流

```cpp
template<class Cap, class Cost>
struct MCMF_Graph {
  public:
    struct Edge {
      int u, v;
      Cap cap, flow;
      Cost cost;
    };
    
    constexpr MCMF_Graph(): n(0) {}
    explicit constexpr MCMF_Graph(int _n): n(_n) {}
    
    constexpr int Add_Edge(int u, int v, Cap cap, Cost cost) {
      assert(0 <= u && u < n);
      assert(0 <= v && v < n);
      assert(0 <= cap);
      assert(0 <= cost);
      
      const int m = static_cast<int>(edge.size());
      edge.emplace_back(u, v, cap, 0, cost);
      return m;
    }
    constexpr Edge Get_Edge(int i) {
      const int m = static_cast<int>(edge.size());
      assert(0 <= i && i < m);
      return edge[i];
    }
    constexpr vector<Edge> Get_All_Edges() {
      return edge;
    }
    constexpr pair<Cap, Cost> MCMF(int S, int T) {
      return MCMF(S, T, numeric_limits<Cap>::max());
    }
    constexpr pair<Cap, Cost> MCMF(int S, int T, Cap lim) {
      return Slope(S, T, lim).back();
    }
    constexpr vector<pair<Cap, Cost>> Slope(int S, int T) {
      return Slope(S, T, numeric_limits<Cap>::max());
    }
    constexpr vector<pair<Cap, Cost>> Slope(int S, int T, Cap lim) {
      assert(0 <= S && S < n);
      assert(0 <= T && T < n);
      assert(S != T);
      
      const int m = static_cast<int>(edge.size());
      vector<int> idx(m), rev(m), deg(n);
      vector<pair<int, _Edge>> elist(2 * m);
      
      for (int i = 0; i < m; i ++) {
        auto e = edge[i];
        idx[i] = deg[e.u] ++;
        rev[i] = deg[e.v] ++;
        elist[i * 2] = {e.u, {e.v, -1, e.cap - e.flow, e.cost}};
        elist[i * 2 + 1] = {e.v, {e.u, -1, e.flow, - e.cost}};
      }
      auto G = Graph(n, elist);
      for (int i = 0; i < m; i ++) {
        auto e = edge[i];
        idx[i] += G.head[e.u];
        rev[i] += G.head[e.v];
        G.elist[idx[i]].rev = rev[i];
        G.elist[rev[i]].rev = idx[i];
      }
      
      auto ans = Slope(G, S, T, lim);
      
      for (int i = 0; i < m; i ++) {
        auto e = G.elist[idx[i]];
        edge[i].flow = edge[i].cap - e.cap;
      }
      
      return ans;
    }
  private:
    int n;
    vector<Edge> edge;
    
    struct _Edge {
      int to, rev;
      Cap cap;
      Cost cost;
    };
    
    struct Graph {
      vector<int> head;
      vector<_Edge> elist;
      
      explicit constexpr Graph(int n,
        const vector<pair<int, _Edge>>& edges):
        head(n + 1), elist(edges.size()) {
          for (auto e : edges) {
            ++ head[e.first + 1];
          }
          for (int i = 1; i <= n; i ++) {
            head[i] += head[i - 1];
          }
          auto cnt = head;
          for (auto e : edges) {
            elist[cnt[e.first] ++] = e.second;
          }
        }
    };
    
    constexpr vector<pair<Cap, Cost>> Slope(Graph& G, int S, int T, Cap lim) {
      struct Element {
        Cost cost;
        int to;
        
        constexpr bool operator < (const Element& rhs) const {
          return cost > rhs.cost;
        }
      };
      
      const Cost inf = numeric_limits<Cost>::max();
      vector<Cost> dual(n), dist(n);
      vector<int> pre(n), Qmn;
      vector<bool> vis(n);
      vector<Element> Q;
      
      auto Dijkstra = [&]() {
        fill(dist.begin(), dist.end(), inf);
        fill(vis.begin(), vis.end(), false);
        Qmn.clear(), Q.clear();
        
        size_t tl = 0;
        
        dist[S] = 0;
        Qmn.emplace_back(S);
        while (Qmn.size() || Q.size()) {
          int u;
          if (Qmn.size()) {
            u = Qmn.back();
            Qmn.pop_back();
          } else {
            while (tl < Q.size()) {
              push_heap(Q.begin(), Q.begin() + (++ tl));
            }
            u = Q.front().to;
            pop_heap(Q.begin(), Q.end());
            Q.pop_back(), -- tl;
          }
          
          if (vis[u]) {
            continue;
          }
          vis[u] = true;
          if (u == T) {
            break;
          }
          
          Cost dual_u = dual[u], dist_u = dist[u];
          for (int i = G.head[u]; i < G.head[u + 1]; i ++) {
            auto e = G.elist[i];
            if (!e.cap) {
              continue;
            }
            Cost cost = e.cost - dual[e.to] + dual_u;
            if (dist_u + cost < dist[e.to]) {
              Cost dist_v = dist_u + cost;
              dist[e.to] = dist_v, pre[e.to] = e.rev;
              if (dist_v == dist_u) {
                Qmn.emplace_back(e.to);
              } else {
                Q.emplace_back(dist_v, e.to);
              }
            }
          }
        }
        
        if (!vis[T]) {
          return false;
        }
        
        for (int u = 0; u < n; u ++) {
          if (!vis[u]) {
            continue;
          }
          dual[u] -= dist[T] - dist[u];
        }
        return true;
      };
      
      Cap flow = 0;
      Cost cost = 0, pre_slope = -1;
      vector<pair<Cap, Cost>> slope = {{Cap(0), Cost(0)}};
      while (flow < lim && Dijkstra()) {
        Cap c = lim - flow;
        Cost d = - dual[S];
        for (int u = T; u != S; u = G.elist[pre[u]].to) {
          c = min(c, G.elist[G.elist[pre[u]].rev].cap);
        }
        for (int u = T; u != S; u = G.elist[pre[u]].to) {
          auto &e = G.elist[pre[u]];
          e.cap += c, G.elist[e.rev].cap -= c;
        }
        flow += c, cost += c * d;
        if (pre_slope == d) {
          slope.back() = {flow, cost};
        } else {
          slope.emplace_back(flow, cost);
        }
        pre_slope = d;
      }
      return slope;
    }
};
```