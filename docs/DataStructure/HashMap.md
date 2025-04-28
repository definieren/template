# 自定义 Hash 函数

```cpp
struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};
```

# 哈希表

```cpp
template<class Key, class T, class Hash = custom_hash, int P = static_cast<int>(1e5 + 7)>
struct Hash_Map {
	private:
		static constexpr int Mod = P;
		static constexpr Hash H{};
		
		struct Info {
			Key key;
			T val;
			int nxt;
		};
		
		int head[Mod], id, vis[Mod];
		vector<Info> D;
	
	public:
		constexpr Hash_Map(): id(1), D{} {}
		template<class _InputIterator, class = _RequireInputIter<_InputIterator>>
		explicit constexpr Hash_Map(_InputIterator __first, _InputIterator __last) {
			D.clear(), id = 1;
			for (auto it = __first; it != __last; it ++) {
				(*this)[it -> first] = it -> second;
			}
		}
		explicit constexpr Hash_Map(const initializer_list<pair<Key, T>> &a) {
			id = 1, D.clear();
			for (auto i : a) {
				dpi(i.fir, i.sec);
				(*this)[i.fir] = i.sec;
			}
		}
		
		constexpr bool Empty() {
			return D.emtpy();
		}
		constexpr size_t Size() {
			return D.size();
		}
		constexpr void Clear() {
			D.clear(), id ++;
		}
		constexpr size_t Count(Key key) {
			int hsh = H(key) % Mod;
			if (vis[hsh] != id) {
				vis[hsh] = id, head[hsh] = -1;
			}
			for (int i = head[hsh]; ~i; i = D[i].nxt) {
				if (D[i].key == key) {
					return 1;
				}
			}
			return 0;
		}
		constexpr vector<pair<Key, T>> Data() {
			const int n = static_cast<int>(D.size());
			vector<pair<Key, T>> data(n);
			for (int i = 0; i < n; i ++) {
				data[i] = {D[i].key, D[i].val};
			}
			return data;
		}
		constexpr void Insert(Key key, T val) {
			(*this)[key] = val;
			return;
		}
		constexpr T operator [] (const Key &key) const {
			int hsh = H(key) % Mod;
			if (vis[hsh] != id) {
				vis[hsh] = id, head[hsh] = -1;
			}
			for (int i = head[hsh]; ~i; i = D[i].nxt) {
				if (D[i].key == key) {
					return D[i].val;
				}
			}
			D.emplace_back(key, T{}, head[key]);
			head[key] = static_cast<int>(D.size()) - 1;
			return D.back().val;
		}
		constexpr T& operator [] (const Key &key) {
			int hsh = H(key) % Mod;
			if (vis[hsh] != id) {
				vis[hsh] = id, head[hsh] = -1;
			}
			for (int i = head[hsh]; ~i; i = D[i].nxt) {
				if (D[i].key == key) {
					return D[i].val;
				}
			}
			D.emplace_back(key, T{}, head[hsh]);
			head[hsh] = static_cast<int>(D.size()) - 1;
			return D.back().val;
		}
};
```