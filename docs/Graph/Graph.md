```cpp
template<class T>
struct Graph {
	private:
		vector<pair<int, T>> edges;
	
	public:
		int n;
		vector<int> head;
		vector<T> elist;
		
		constexpr Graph(): n(0), edges{} {}
		explicit constexpr Graph(int _n): n(_n), edges{} {}
		explicit constexpr Graph(int _n, const vector<pair<int, T>> &es) {
			_n = n, edges = es, Build();
		}
		constexpr void Add_Edge(int u, T w) {
			edges.push_back({u, w});
			return;
		}
		constexpr void Build() {
			head.assign(n + 1, 0), elist.assign(edges.size(), T{});
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
			return;
		}
		template<class F = void(*)(T&)>
		constexpr void Trans(int u, F f) {
			for (int i = head[u]; i < head[u + 1]; i ++) {
				f(elist[i]);
			}
			return;
		}
		constexpr bool Empty(int u) {
			return head[u] == head[u + 1];
		}
};
```