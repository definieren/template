用于离线快速存图。

`graph<T, directed>`：`T` 为边权类型，`directed` 表示是否有向。

可以通过 `addEdge(u, v, cost, id)` 加边，`build()` 离线建图。

`operator []` 提供了遍历一个点出边的方式。

`degArray(), IOdegArray(), deg(), indeg(), outdeg()` 可以求出每个点的度数、入度、出度。

```cpp
namespace graphBase {
template<class T>
struct edge {
  int u, v;
  T cost;
  int id;
  
  edge rev() const {
    return {v, u, cost, id};
  }
};
template<class T = int, bool directed = false>
struct graph {
public:
  static constexpr bool is_directed = directed;
  int n, m;
  std::vector<edge<T>> edges;
  std::vector<int> head;
  std::vector<edge<T>> csr_edges;
  
private:
  std::vector<int> deg_, indeg_, outdeg_;
  bool prepared;
  
  struct edgeList {
    using iter = std::vector<edge<T>>::iterator;
  private:
    iter begin_, end_;
  public:
    edgeList(iter begin, iter end): begin_(begin), end_(end) {}
    iter begin() const {
      return begin_;
    }
    iter end() const {
      return end_;
    }
    auto size() const {
      return std::distance(begin_, end_);
    }
    bool empty() const {
      return !size();
    }
    edge<T>& operator [] (int i) const {
      return begin_[i];
    }
  };
  struct edgeListConst {
    using iter = std::vector<edge<T>>::const_iterator;
  private:
    iter begin_, end_;
  public:
    edgeListConst(iter begin, iter end): begin_(begin), end_(end) {}
    iter begin() const {
      return begin_;
    }
    iter end() const {
      return end_;
    }
    auto size() const {
      return std::distance(begin_, end_);
    }
    bool empty() const {
      return !size();
    }
    const edge<T>& operator [] (int i) const {
      return begin_[i];
    }
  };
  
public:
  bool is_prepared() const {
    return prepared;
  }
  
  graph(): n(0), m(0), deg_{}, indeg_{}, outdeg_{}, prepared(false) {}
  graph(int n): n(n), m(0), deg_{}, indeg_{}, outdeg_{}, prepared(false) {}
  graph(int n_, const std::vector<edge<T>>& es, const int offset = 1):
    n(n_), m(es.size()), edges(es),
    deg_{}, indeg_{}, outdeg_{}, prepared(false) {
    for (auto& e : edges) {
      e.u -= offset, e.v -= offset;
    }
    build();
  }
  
  void addEdge(int u, int v, T cost = 1, int id = -1) {
    assert(!prepared);
    assert(0 <= u && u < n);
    assert(0 <= v && v < n);
    if (id == -1) {
      id = m;
    }
    edges.emplace_back(u, v, cost, id);
    ++ m;
  }
  
  void build() {
    assert(!prepared);
    prepared = true;
    head.assign(n + 1, 0);
    for (auto&& e : edges) {
      ++ head[e.u + 1];
      if constexpr (!is_directed) {
        ++ head[e.v + 1];
      }
    }
    for (int u = 0; u < n; u ++) {
      head[u + 1] += head[u];
    }
    auto ptr = head;
    csr_edges.resize(head.back() + 1);
    for (auto&& e : edges) {
      csr_edges[ptr[e.u] ++] = e;
      if constexpr (!is_directed) {
        csr_edges[ptr[e.v] ++] = e.rev();
      }
    }
  }
  
  edgeList operator [] (int u) {
    assert(prepared);
    assert(0 <= u && u < n);
    return {csr_edges.begin() + head[u], csr_edges.begin() + head[u + 1]};
  }
  edgeListConst operator [] (int u) const {
    assert(prepared);
    assert(0 <= u && u < n);
    return {csr_edges.begin() + head[u], csr_edges.begin() + head[u + 1]};
  }
  
private:
  void degArray_() {
    assert(prepared);
    assert(deg_.empty());
    deg_.resize(n);
    for (auto&& e : edges) {
      ++ deg_[e.u], ++ deg_[e.v];
    }
  }
  void IOdegArray_() {
    assert(prepared);
    assert(indeg_.empty());
    assert(outdeg_.empty());
    indeg_.resize(n);
    outdeg_.resize(n);
    for (auto&& e : edges) {
      ++ indeg_[e.v];
      ++ outdeg_[e.u];
    }
  }
  
public:
  std::vector<int> degArray() {
    if (deg_.empty()) {
      degArray_();
    }
    return deg_;
  }
  std::pair<std::vector<int>, std::vector<int>> IOdegArray() {
    if (indeg_.empty()) {
      IOdegArray_();
    }
    return {indeg_, outdeg_};
  }
  int deg(int u) {
    assert(0 <= u && u < n);
    if (deg_.empty()) {
      degArray_();
    }
    return deg_[u];
  }
  int indeg(int u) {
    assert(0 <= u && u < n);
    if (indeg_.empty()) {
      IOdegArray_();
    }
    return indeg_[u];
  }
  int outdeg(int u) {
    assert(0 <= u && u < n);
    if (outdeg_.empty()) {
      IOdegArray_();
    }
    return outdeg_[u];
  }
};
} using namespace graphBase;
```