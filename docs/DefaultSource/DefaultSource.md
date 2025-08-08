简略版：

```cpp
#include <bits/stdc++.h>

using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using i128 = __int128;
using u128 = unsigned __int128;

int main() {
#ifdef LOCAL
	freopen("!in.in", "r", stdin);
	freopen("!out.out", "w", stdout);
#endif
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	std::cout.tie(nullptr);



	return 0;
}
```

长版：

```cpp
#include <bits/stdc++.h>

#define fir first
#define sec second
#define mkp make_pair
#define mkt make_tuple
#ifdef LOCAL
#define dbg(x) cerr << "In Line " << __LINE__ << " the " << #x << " = " << x << '\n'
#define dpi(x, y) cerr << "In Line " << __LINE__ << " the " << #x << " = " << x << " ; " << "the " << #y << " = " << y << '\n'
#define dbgf(fmt, args...) fprintf(stderr, fmt, ##args)
#else
#define dbg(x) void()
#define dpi(x, y) void()
#define dbgf(fmt, args...) void()
#endif

using namespace std;

using i64 = long long;
using u64 = unsigned long long;
using u32 = unsigned int;
using ldb = long double;
using i128 = __int128_t;
using u128 = __uint128_t;
using pii = pair<int, int>;
using pil = pair<int, i64>;
using pli = pair<i64, int>;
using vi = vector<int>;
using vpii = vector<pii>;

namespace {
bool Mbe;
constexpr int MOD = 998244353;
template<typename T> T Norm(T a, T p = MOD) { return (a % p + p) % p; }
template<typename T> bool cmax(T &a, T b) { return a < b ? a = b, true : false; }
template<typename T> bool cmin(T &a, T b) { return a > b ? a = b, true : false; }
template<typename T> T DivFloor(T a, T b) { return a >= 0 ? a / b : (a - b + 1) / b; }
template<typename T> T DivCeil(T a, T b) { return a >= 0 ? (a + b - 1) / b : a / b; }

namespace FastIO {
	constexpr int LEN = 1 << 20;
	char in[LEN + 1], out[LEN + 1];
	char *pin = in, *pout = out, *ein = in, *eout = out + LEN;

	char gc() { return pin == ein && (ein = (pin = in) + fread(in, 1, LEN, stdin), ein == in) ? EOF : *pin ++; }
	void pc(char c) { pout == eout && (fwrite(out, 1, LEN, stdout), pout = out); (*pout ++) = c; return; }
	struct Flush { ~Flush() { fwrite(out, 1, pout - out, stdout); pout = out; return; } } _flush;

	template<typename T> T Read() {
		T x = 0; int f = 1; char ch = gc();
		while (ch < '0' || ch > '9') f = (ch == '-' ? (~f + 1) : f), ch = gc();
		while (ch >= '0' && ch <= '9') x = (x << 1) + (x << 3) + (ch ^ 48), ch = gc();
		return x * f;
	}
	void Read(char *s) {
		char ch = gc();
		while (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') ch = gc();
		while ((ch != EOF) && !(ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t')) *s = ch, s ++, ch = gc();
		*s = '\0'; return;
	}
	template<typename T> void Read(T &x) { x = Read<T>(); return; }
	template<typename T, typename ...Args>
	void Read(T &x, Args &...args) { Read(x), Read(args...); return; }
	template<typename T> void Write(T x) {
		static char stk[40]; int tp = 0;
		if (x < 0) pc('-'), x = ~x + 1;
		do stk[tp++] = x % 10 + 48, x /= 10; while (x);
		while (tp --) pc(stk[tp]);
		return;
	}
	void Write(char ch) { pc(ch); return; }
	void Write(const char *s) {
		while (*s != '\0') pc(*s), s ++;
		return;
	}
	void Puts(const char *s) {
		Write(s), pc('\n'); return;
	}
	template<typename T, typename ...Args>
	void Write(T x, Args ...args) { Write(x), Write(args...); return; }
} using namespace FastIO;

void slv() {
	
	return;
}
void clr() {
	
	return;
}
bool Med;
}

int main() {
#ifdef LOCAL
	freopen("!in.in", "r", stdin);
	freopen("!out.out", "w", stdout);
	fprintf(stderr, "%.3lf Mb\n", fabs((&Mbe - &Med) / 1048576.0));
#endif
	int T = 1;
//	int T = Read<int>();
	while (T --) slv(), clr();
#ifdef LOCAL
	fprintf(stderr, "%d ms\n", (int)clock());
#endif
	return 0;
}
```