```cpp
template<class Info, class Tag>
struct Segment_Tree {
	private:
		int n;
		vector<int> ps;
		struct Node { Info dt; Tag tg; };
		vector<Node> sgt;
		vector<Info> ve; vector<int> _vec;
		vector<tuple<int, int, int>> vec;
		
		constexpr void Push_Up(int p) { sgt[p].dt = sgt[p << 1].dt + sgt[p << 1 | 1].dt; return; }
		constexpr void Push_Tag(int p, Tag tg) { sgt[p].dt = tg.Apply(sgt[p].dt), sgt[p].tg = tg * sgt[p].tg; return; }
		constexpr void Push_Down(int p) { if (sgt[p].tg) Push_Tag(p << 1, sgt[p].tg), Push_Tag(p << 1 | 1, sgt[p].tg), sgt[p].tg = Tag(); return; }
		constexpr void Build(int p, int L, int R) {
			sgt[p].tg = Tag();
			if (L == R) { ps[L] = p; sgt[p].dt = ve[L]; return; }
			int Mid = (L + R) >> 1; Build(p << 1, L, Mid);
			Build(p << 1 | 1, Mid + 1, R); Push_Up(p); return;
		}
		constexpr void Extract(int l, int r, int p, int L, int R) {
			if (l <= L && R <= r) { vec.emplace_back(p, L, R); return; }
			Push_Down(p); int Mid = (L + R) >> 1;
			if (l <= Mid) Extract(l, r, p << 1, L, Mid);
			if (Mid < r) Extract(l, r, p << 1 | 1, Mid + 1, R);
			_vec.emplace_back(p); return;
		}
	public:
		constexpr Segment_Tree(int _n) {
			Build(_n);
		}
		template<class T>
		constexpr Segment_Tree(const vector<T> &a) {
			Build(a);
		}
		template<class _InputIterator, class = _RequireInputIter<_InputIterator>>
		constexpr Segment_Tree(_InputIterator __first, _InputIterator __last) {
			Build(__first, __last);
		}
		template<class F>
		constexpr Segment_Tree(int _n, F f) {
			Build<F>(_n, f);
		}
		constexpr void Build(int _n) {
			ve.assign(n = _n, Info{});
			sgt.assign(n << 2, Node{});
			ps.assign(n, 0);
			Build(1, 0, n - 1); return;
		}
		template<class T> constexpr void Build(const vector<T> &a) {
			ve.resize(n = (int)a.size());
			sgt.assign(n << 2, Node{});
			ps.assign(n, 0);
			for (int i = 0; i < n; i ++) ve[i] = a[i];
			Build(1, 0, n - 1); return;
		}	
		template<class _InputIterator, class = _RequireInputIter<_InputIterator>>
		constexpr void Build(_InputIterator __first, _InputIterator __last) {
			ve.resize(n = __last - __first);
			sgt.assign(n << 2, Node{});
			ps.assign(n, 0);
			for (auto it = __first; it != __last; it ++) ve[it - __first] = *it;
			Build(1, 0, n - 1); return;
		}
		template<class F> constexpr void Build(int _n, F f) {
			ve.resize(n = _n);
			sgt.assign(n << 2, Node{});
			ps.assign(n, 0);
			for (int i = 0; i < n; i ++) ve[i] = f(i);
			Build(1, 0, n - 1); return;
		}
		constexpr void Update(int l, int r, Tag k) {
			_vec.clear(), vec.clear();
			Extract(l, r, 1, 0, n - 1);
			for (auto [p, L, R] : vec) Push_Tag(p, k);
			for (auto p : _vec) Push_Up(p);
			return;
		}
		constexpr void Update(int pos, Tag k) {
			_vec.clear(); pos = ps[pos];
			while (pos) _vec.emplace_back(pos), pos >>= 1;
			for (int i = (int)_vec.size() - 1; i; i --) Push_Down(_vec[i]);
			Push_Tag(_vec.front(), k);
			for (int i = 1; i < (int)_vec.size(); i ++) Push_Up(_vec[i]);
			return;
		}
		constexpr void Change(int pos, Info k) {
			_vec.clear(); pos = ps[pos];
			while (pos) _vec.emplace_back(pos), pos >>= 1;
			for (int i = (int)_vec.size() - 1; i; i --) Push_Down(_vec[i]);
			sgt[_vec.front()].dt = k;
			for (int i = 1; i < (int)_vec.size(); i ++) Push_Up(_vec[i]);
			return;
		}
		constexpr Info Query() { return sgt[1].dt; }
		constexpr Info Query(int l, int r) {
			_vec.clear(), vec.clear();
			Extract(l, r, 1, 0, n - 1);
			Info ans = sgt[get<0>(vec.front())].dt;
			for (int i = 1; i < (int)vec.size(); i ++)
				ans = ans + sgt[get<0>(vec[i])].dt;
			return ans;
		}
		constexpr Info Query(int pos) {
			_vec.clear(); pos = ps[pos];
			while (pos) _vec.emplace_back(pos), pos >>= 1;
			for (int i = (int)_vec.size() - 1; i; i --) Push_Down(_vec[i]);
			return sgt[_vec.front()].dt;
		}
		template<class Compare>
		constexpr pair<int, Info> Find_Nxt(int l, int r, Info k, Compare comp) {
			_vec.clear(), vec.clear();
			Extract(l, r, 1, 0, n - 1);
			for (auto [p, L, R] : vec) {
				if (!comp(k, sgt[p].dt)) continue;
				while (L < R) {
					int Mid = (L + R) >> 1; Push_Down(p);
					if (comp(k, sgt[p << 1].dt)) R = Mid, p <<= 1;
					else L = Mid + 1, p = p << 1 | 1;
				}
				return {L, sgt[p].dt};
			}
			return {r + 1, Info()};
		}
		template<class Compare>
		constexpr pair<int, Info> Find_Pre(int l, int r, Info k, Compare comp) {
			_vec.clear(), vec.clear();
			Extract(l, r, 1, 0, n - 1);
			reverse(vec.begin(), vec.end());
			for (auto [p, L, R] : vec) {
				if (!comp(k, sgt[p].dt)) continue;
				while (L < R) {
					int Mid = (L + R) >> 1; Push_Down(p);
					if (!comp(k, sgt[p << 1 | 1].dt)) R = Mid, p <<= 1;
					else L = Mid + 1, p = p << 1 | 1;
				}
				return {L, sgt[p].dt};
			}
			return {l - 1, Info()};
		}
};
struct Info {
	
	Info() {}
	friend Info operator + (Info x, Info y) {
		
	}
};
struct Tag {
	
	Tag() {}
	operator bool() {
		
	}
	friend Tag operator * (Tag x, Tag y) {
		
	}
	Info Apply(Info x) {
		
	}
};
```