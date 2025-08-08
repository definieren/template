并查集，使用路径压缩加按秩合并，单次均摊 $O(\alpha(n))$。

```cpp
struct DisjointSet {
private:
	int n;
	std::vector<int> fa;
	
	int _find(int u) {
		if (fa[u] < 0) {
			return u;
		}
		return fa[u] = _find(fa[u]);
	}
public:
	DisjointSet(): n{}, fa{} {}
	DisjointSet(int n_): n(n_), fa(n, -1) {}
	
	int find(int u) {
		assert(0 <= u && u < n);
		return _find(u);
	}
	int merge(int u, int v) {
		assert(0 <= u && u < n);
		assert(0 <= v && v < n);
		u = _find(u), v = _find(v);
		if (u == v) {
			return u;
		}
		if (-fa[u] < -fa[v]) {
			std::swap(u, v);
		}
		fa[u] += fa[v], fa[v] = u;
		return u;
	}
	int size(int u) {
		assert(0 <= u && u < n);
		return -fa[_find(u)];
	}
	bool same(int u, int v) {
		assert(0 <= u && u < n);
		assert(0 <= v && v < n);
		return _find(u) == _find(v);
	}
	std::vector<std::vector<int>> groups() {
		std::vector<std::vector<int>> ans(n);
		for (int u = 0; u < n; u ++) {
			if (u == _find(u)) {
				ans[u].reserve(-fa[u]);
			}
		}
		for (int u = 0; u < n; u ++) {
			ans[_find(u)].push_back(u);
		}
		ans.erase(std::remove_if(ans.begin(), ans.end(),
			[&](const std::vector<int>& v) { return v.empty(); }), ans.end());
		return ans;
	}
};
```