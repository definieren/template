# 双极定向

对于图 $G = (V, E)$，以及 $s, t \in V$，以下四个命题等价：

- 图 $G' = (V, E \cup \{ (s, t) \})$ 点双连通。
- 圆方树上所有方点成链，并且 $s \rightsquigarrow t$ 是圆方树的直径。
- 存在 DAG $G'$ 以 $G$ 为基图，且 $s$ 是唯一入度为 $0$ 的点，$t$ 是唯一出度为 $0$ 的点。
- 存在 $p \in \mathfrak S_n$，$p_0 = s, p_{n - 1} = t$，任意前后缀导出子图连通。

`st_Numbering(G, s, t)`：返回图 $G$ 中满足 $p_0 = s, p_{n - 1} = t$ 且任意前后缀导出子图连通的排列 $p$。

`vector` 存图版本：

```cpp
std::vector<int> st_Numbering(const std::vector<std::vector<int>>& adj, const int& s, const int& t) {
	const int n = adj.size();
	if (n == 1) {
		return {0};
	}
	if (s == t) {
		return {};
	}
	std::vector<int> dfn(n, -1), fa(n, -1), low(n, -1), ord(n);
	{
		std::vector<u32> head(n);
		int u = t, v, dfc = 2;
		bool reachable = false;
		ord[0] = s, ord[1] = t;
		low[s] = dfn[s] = 0;
		low[t] = dfn[t] = 1;
		fa[t] = s;
		while (u != s) {
			if (head[u] == adj[u].size()) {
				v = u, u = fa[v];
				if (~u) {
					low[u] = std::min(low[u], low[v]);
				}
				continue;
			}
			v = adj[u][head[u] ++];
			if (v == s) {
				reachable = true;
			}
			if (dfn[v] == -1) {
				ord[dfc] = v;
				low[v] = dfn[v] = dfc ++;
				fa[v] = u, u = v;
				continue;
			}
			low[u] = std::min(low[u], dfn[v]);
		}
		
		if (!reachable || dfc < n) {
			return {};
		}
	};
	
	std::vector<bool> sgn(n);
	std::vector<int> L(n), R;
	std::iota(L.begin(), L.end(), 0), R = L;
	for (int i = 1; i < n; i ++) {
		const int u = ord[i];
		if (u != t && low[u] == dfn[fa[u]]) {
			return {};
		}
		if (sgn[ord[low[u]]]) {
			const int v = L[fa[u]];
			L[u] = v, R[v] = u;
			R[u] = fa[u], L[fa[u]] = u;
			sgn[fa[u]] = false;
		} else {
			const int v = R[fa[u]];
			R[u] = v, L[v] = u;
			L[u] = fa[u], R[fa[u]] = u;
			sgn[fa[u]] = true;
		}
	}
	
	std::vector<int> p(n);
	ord[0] = s;
	for (int i = 1; i < n; i ++) {
		ord[i] = R[ord[i - 1]];
	}
	for (int i = 0; i < n; i ++) {
		p[ord[i]] = i;
	}
	
	return p;
}
```

`graph` 存图版本（需前置 `GraphBase`）：

```cpp
template<class Graph>
std::vector<int> st_Numbering(const Graph& G, const int& s, const int& t) {
	assert(!G.is_directed);
	assert(G.is_prepared());
	const int n = G.n;
	if (n == 1) {
		return {0};
	}
	if (s == t) {
		return {};
	}
	std::vector<int> dfn(n, -1), fa(n, -1), low(n, -1), ord(n);
	{
		auto head = G.head;
		int u = t, v, dfc = 2;
		bool reachable = false;
		ord[0] = s, ord[1] = t;
		low[s] = dfn[s] = 0;
		low[t] = dfn[t] = 1;
		fa[t] = s;
		while (u != s) {
			if (head[u] == G.head[u + 1]) {
				v = u, u = fa[v];
				if (~u) {
					low[u] = std::min(low[u], low[v]);
				}
				continue;
			}
			v = G.csr_edges[head[u] ++].v;
			if (v == s) {
				reachable = true;
			}
			if (dfn[v] == -1) {
				ord[dfc] = v;
				low[v] = dfn[v] = dfc ++;
				fa[v] = u, u = v;
				continue;
			}
			low[u] = std::min(low[u], dfn[v]);
		}
		
		if (!reachable || dfc < n) {
			return {};
		}
	};
	
	std::vector<bool> sgn(n);
	std::vector<int> L(n), R;
	std::iota(L.begin(), L.end(), 0), R = L;
	for (int i = 1; i < n; i ++) {
		const int u = ord[i];
		if (u != t && low[u] == dfn[fa[u]]) {
			return {};
		}
		if (sgn[ord[low[u]]]) {
			const int v = L[fa[u]];
			L[u] = v, R[v] = u;
			R[u] = fa[u], L[fa[u]] = u;
			sgn[fa[u]] = false;
		} else {
			const int v = R[fa[u]];
			R[u] = v, L[v] = u;
			L[u] = fa[u], R[fa[u]] = u;
			sgn[fa[u]] = true;
		}
	}
	
	std::vector<int> p(n);
	ord[0] = s;
	for (int i = 1; i < n; i ++) {
		ord[i] = R[ord[i - 1]];
	}
	for (int i = 0; i < n; i ++) {
		p[ord[i]] = i;
	}
	
	return p;
}
```