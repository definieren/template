```cpp
struct BigNum {
	private:
		static constexpr ll lim = 1e8, mod = 1e16, lm = 10,
			Mod = 167772161, g = 3, inv2 = (Mod + 1) >> 1;
		static constexpr int pw10[] = {1, 10, 100, 1000,
			10000, 100000, 1000000, 10000000, 100000000};
		static vector<char> str;
		static vector<int> rev;
		static vector<ll> w;
		
		i64 Mul(i64 x, i64 y) { return x * y % Mod; }
		i64 Pow(i64 a, i64 b) {
			i64 ans = 1;
			for (; b; b >>= 1, a = Mul(a, a))
				if (b & 1) ans = Mul(ans, a);
			return ans;
		}
		void DFT(vector<i64> &a) {
			const int n = (int)a.size();
			if ((int)rev.size() != n) { 
				int k = __lg(n) - 1; rev.resize(n);
				for (int i = 0; i < n; i ++)
					rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << k);
			}
			if ((int)w.size() < n) {
				int k = __lg(w.size()); w.resize(n);
				while (1 << k < n) {
					i64 e = Pow(g, (Mod - 1) >> (k + 1));
					for (int i = 1 << (k - 1); i < 1 << k; i ++)
						w[i << 1] = w[i], w[i << 1 | 1] = Mul(w[i], e);
					k ++;
				}
			}
			for (int i = 0; i < n; i ++)
				if (i < rev[i]) swap(a[i], a[rev[i]]);
			for (int k = 1; k < n; k <<= 1)
				for (int i = 0; i < n; i += k << 1)
					for (int j = 0; j < k; j ++) {
						i64 u = a[i | j], v = Mul(a[i | j | k], w[j | k]);
						a[i | j] = u + v >= Mod ? u + v - Mod : u + v;
						a[i | j | k] = u - v < 0 ? u - v + Mod : u - v;
					}
			return;
		}
		void IDFT(vector<i64> &a) {
			const int n = (int)a.size();
			reverse(a.begin() + 1, a.end());
			DFT(a); i64 inv = (1 - Mod) / n + Mod;
			for (int i = 0; i < n; i ++) a[i] = Mul(a[i], inv);
			return;
		}
		vector<i64> Convolution(vector<i64> f, vector<i64> g) {
			if (f.size() < g.size()) std::swap(f, g);
			const int N = f.size(), M = g.size();
			int n = 1, m = N + M - 1;
			if (g.size() < 1 << 7) {
				vector<i64> ret(m);
				for (int i = 0; i < N; i ++)
					for (int j = 0; j < M; j ++)
						ret[i + j] += f[i] * g[j];
				return ret;
			}
			while (n < m) n <<= 1;
			f.resize(n), g.resize(n), DFT(f), DFT(g);
			for (int i = 0; i < n; i ++) f[i] = Mul(f[i], g[i]);
			IDFT(f), f.resize(m);
			return f;
		}
	public:
		bool sign;
		vector<i64> num;
		
		BigNum(): sign(false), num(1, 0) { return; }
		BigNum(i64 x) {
			if (x < 0) sign = true, x = - x; else sign = false;
			do num.emplace_back(x % lim), x /= lim; while (x);
			return;
		}
		void Read() {
			char ch = getchar();
			sign = false, num.clear(), str.clear();
			while (ch < '0' || ch > '9')
				sign ^= ch == '-', ch = getchar();
			while ('0' <= ch && ch <= '9')
				str.emplace_back(ch), ch = getchar();
			const int n = (int)str.size(); num.resize((n + 7) / 8);
			reverse(str.begin(), str.end());
			for (int i = 0, j = -1; i < n; i ++) {
				if (i % 8 == 0) j ++;
				num[j] += pw10[i - j * 8] * (str[i] - '0');
			}
			return;
		}
		void Write() {
			if (!(*this)) { putchar('0'); return; }
			if (sign) putchar('-');
			const int n = num.size();
			::Write(num[n - 1]);
			for (int i = n - 2; ~i; i --)
				for (int j = 7; ~j; j --)
					putchar(num[i] / pw10[j] % 10 + '0');
			return;
		}
		vector<i64> To_Base_1() const {
			const int n = num.size();
			vector<i64> ret(n << 3);
			for (int i = 0; i < n; i ++) {
				i64 x = num[i];
				for (int j = 0; j < 8; j ++)
					ret[i * 8 + j] = x % lm, x /= lm;
			}
			return ret;
		}
		void To_Base_8(vector<i64> v) {
			const int n = v.size();
			num.assign((n + 7) >> 3, 0);
			for (int i = 0, j = -1; i < n; i ++) {
				if (!(i & 7)) j ++;
				num[j] += v[i] * pw10[i - j * 8];
			}
			while (num.size() > 1 && !num.back())
				num.pop_back();
			return;
		}
		explicit operator bool() const {
			return num.size() != 1 || num[0];
		}
		int get(const int &i) const {
			int x = i / 8, y = i % 8;
			return num[x] / pw10[y] % 10;
		}
		BigNum Abs() const {
			BigNum tmp = *this; tmp.sign = false;
			return tmp;
		}
		BigNum Shift(int k) {
			BigNum tmp = *this;
			tmp.num.insert(tmp.num.begin(), k, 0);
			return tmp;
		}
		friend BigNum operator - (BigNum x) {
			return x.sign ^= 1, x;
		}
		BigNum& operator += (const BigNum& oth) {
			if (sign == oth.sign) {
				const int m = oth.num.size(), n = max((int)num.size(), m);
				num.resize(n, 0);
				for (int i = 0; i < m; i ++) num[i] += oth.num[i];
				for (int i = 0; i + 1 < n; i ++)
					num[i + 1] += num[i] >= lim ? 1 : 0,
					num[i] = num[i] >= lim ? num[i] - lim : num[i];
				if (num[n - 1] >= lim)
					num.emplace_back(1), num[n - 1] -= lim;
				return *this;
			}
			if (sign) return (*this) = oth - (- (*this));
			return (*this) -= - oth;
		}
		BigNum& operator -= (const BigNum& oth) {
			if (oth.sign == sign) {
				if ((*this) == oth) return (*this) = BigNum(0);
				if ((*this).Abs() < oth.Abs()) return (*this) = - (oth - (*this));
				const int n = num.size(), m = oth.num.size();
				for (int i = 0; i < m; i ++) num[i] -= oth.num[i];
				for (int i = 0; i + 1 < n; i ++)
					if (num[i] < 0) num[i] += lim, num[i + 1] --;
				while (num.size() > 1 && !num.back())
					num.pop_back();
				return *this;
			}
			if (sign) return (*this) = - ((- (*this)) + oth);
			return (*this) += - oth;
		}
		BigNum& operator *= (const BigNum& oth) {
			if (!(*this)) return *this;
			if (!oth) return (*this) = oth;
			const int n = num.size(), m = oth.num.size();
			if ((m < 1 << 7) || (n < 1 << 7)) {
				BigNum ans; ans.sign = sign ^ oth.sign;
				ans.num.resize(n + m, 0);
				for (int i = 0; i < n; i ++)
					for (int j = 0; j < m; j ++) {
						ans.num[i + j] += num[i] * oth.num[j];
						if (ans.num[i + j] >= mod)
							ans.num[i + j + 1] += lim, ans.num[i + j] -= mod;
					}
				for (int i = 0; i + 1 < n + m; i ++)
					ans.num[i + 1] += ans.num[i] / lim,
					ans.num[i] %= lim;
				if (!ans.num.back()) ans.num.pop_back();
				return (*this) = ans;
			}
			BigNum ans; ans.sign = sign ^ oth.sign;
			ans.num = Convolution(To_Base_1(), oth.To_Base_1());
			ans.num.resize((n + m) << 3, 0);
			const int k = ans.num.size();
			for (int i = 0; i + 1 < k; i ++)
				ans.num[i + 1] += ans.num[i] / lm,
				ans.num[i] %= lm;
			ans.To_Base_8(ans.num);
			return (*this) = ans;
		}
		BigNum& operator /= (const BigNum& oth) {
			if (!(*this)) return *this;
			const int n = num.size(), m = oth.num.size();
			BigNum ans, x = Abs(), y = oth.Abs();
			ans.num.resize(max(n - m + 1, 1));
			for (int i = max(n - m, - 1); ~i; i --) {
				int l = 0, r = lim - 1;
				while (l < r) {
					int mid = (l + r + 1) >> 1;
					if ((y * mid).Shift(i) <= x) l = mid;
					else r = mid - 1;
				}
				ans.num[i] = l, x -= (y * l).Shift(i);
			}
			while (ans.num.size() > 1 && !ans.num.back())
				ans.num.pop_back();
			ans.sign = sign ^ oth.sign;
			if (ans.sign && ans * oth != (*this)) -- ans;
			if (ans.num.size() == 1 && !ans.num[0]) ans.sign = false;
			return (*this) = ans;
		}
		BigNum& operator %= (const BigNum& oth) {
			return (*this) = (*this) - (*this) / oth * oth;
		}
		BigNum& operator *= (const i64 &oth) {
			if (!oth) return (*this) = BigNum(0);
			for (auto &x : num) x *= oth;
			for (int i = 0; i + 1 < (int)num.size(); i ++)
				num[i + 1] += num[i] / lim, num[i] %= lim;
			while (num.back() > lim) {
				i64 x = num.back() / lim;
				num.back() %= lim, num.emplace_back(x);
			}
			while (num.size() > 1 && !num.back())
				num.pop_back();
			return *this;
		}
		BigNum& operator /= (const i64 &_x) {
			const int n = num.size();
			BigNum tmp = *this; i64 x = _x;
			if (x < 0) sign ^= 1, x = - x;
			for (int i = n - 1; i; i --)
				num[i - 1] += num[i] % x * lim, num[i] /= x;
			num[0] /= x;
			while (num.size() > 1 && !num.back())
				num.pop_back();
			if (sign && (*this) * x != tmp) -- (*this);
			if (num.size() == 1 && !num[0]) sign = false;
			return *this;
		}
		BigNum& operator %= (const i64 &x) {
			return (*this) = (*this) - (*this) / x * x;
		}
		friend BigNum operator + (BigNum x, const BigNum &y) { return x += y; }
		friend BigNum operator - (BigNum x, const BigNum &y) { return x -= y; }
		friend BigNum operator * (BigNum x, const BigNum &y) { return x *= y; }
		friend BigNum operator / (BigNum x, const BigNum &y) { return x /= y; }
		friend BigNum operator % (BigNum x, const BigNum &y) { return x %= y; }
		friend BigNum operator + (BigNum x, const i64 &y) { return x += BigNum(y); }
		friend BigNum operator - (BigNum x, const i64 &y) { return x -= BigNum(y); }
		friend BigNum operator * (BigNum x, const i64 &y) { return x *= y; }
		friend BigNum operator / (BigNum x, const i64 &y) { return x /= y; }
		friend BigNum operator % (BigNum x, const i64 &y) { return x %= y; }
		friend BigNum operator + (i64 x, const BigNum &y) { return BigNum(x) += y; }
		friend BigNum operator - (i64 x, const BigNum &y) { return BigNum(x) -= y; }
		friend BigNum operator * (const i64 &x, BigNum y) { return y *= x; }
		friend BigNum operator / (i64 x, const BigNum &y) { return BigNum(x) /= y; }
		friend BigNum operator % (i64 x, const BigNum &y) { return BigNum(x) %= y; }
		BigNum operator ++ (int) { BigNum tmp = *this; (*this) += 1; return tmp; }
		BigNum operator -- (int) { BigNum tmp = *this; (*this) -= 1; return tmp; }
		friend BigNum& operator ++ (BigNum& x) { return x += 1; }
		friend BigNum& operator -- (BigNum& x) { return x -= 1; }
		friend bool operator == (const BigNum &x, const BigNum &y) {
			const int n = x.num.size(), m = y.num.size();
			if (x.sign != y.sign || n != m) return false;
			for (int i = 0; i < n; i ++)
				if (x.num[i] != y.num[i]) return false;
			return true;
		}
		friend bool operator != (const BigNum &x, const BigNum &y) {
			return !(x == y);
		}
		friend bool operator < (const BigNum &x, const BigNum &y) {
			if (x.sign ^ y.sign) return x.sign;
			if (x.sign && y.sign) return x.Abs() > y.Abs();
			const int n = x.num.size(), m = y.num.size();
			if (n != m) return n < m;
			for (int i = n - 1; ~i; i --)
				if (x.num[i] != y.num[i]) return x.num[i] < y.num[i];
			return false;
		}
		friend bool operator > (const BigNum &x, const BigNum &y) {
			if (x.sign ^ y.sign) return y.sign;
			if (x.sign && y.sign) return x.Abs() < y.Abs();
			const int n = x.num.size(), m = y.num.size();
			if (n != m) return n > m;
			for (int i = n - 1; ~i; i --)
				if (x.num[i] != y.num[i]) return x.num[i] > y.num[i];
			return false;
		}
		friend bool operator <= (const BigNum &x, const BigNum &y) {
			return !(x > y);
		}
		friend bool operator >= (const BigNum &x, const BigNum &y) {
			return !(x < y);
		}
		friend ostream& operator << (ostream& os, const BigNum &x) {
			if (!x) { os << 0; return os; }
			if (x.sign) os << '-';
			const int n = x.num.size();
			os << x.num[n - 1];
			for (int i = n - 2; ~i; i --)
				for (int j = 7; ~j; j --)
					os << (char)(x.num[i] / pw10[j] % 10 + '0');
			return os;
		}
};
vector<char> BigNum::str;
constexpr int BigNum::pw10[];
vector<int> BigNum::rev;
vector<i64> BigNum::w{0, 1};
```