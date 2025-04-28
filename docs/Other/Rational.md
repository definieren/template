```cpp
template<class T> struct Frac {
	T x, y;
	
	Frac(): x(0), y(1) { return; }
	Frac(T _x): x(_x), y(1) { return; }
	Frac(T _x, T _y): x(_x), y(_y) { if (y < 0) x = - x, y = - y; return; }
	explicit operator double() const { return 1. * x / y; }
	explicit operator bool() const { return x; }
	void Reduce() { T g = __gcd(x, y); x /= g, y /= g; return; }
	Frac& operator += (const Frac &oth) { x = x * oth.y + oth.x * y, y *= oth.y; return *this; }
	Frac& operator -= (const Frac &oth) { x = x * oth.y - oth.x * y, y *= oth.y; return *this; }
	Frac& operator *= (const Frac &oth) { x *= oth.x, y *= oth.y; return *this; }
	Frac& operator /= (const Frac &oth) { x *= oth.y, y *= oth.x; if (y < 0) x = - x, y = - y; return *this; }
	friend Frac operator + (Frac x, const Frac &y) { return x += y; }
	friend Frac operator - (Frac x, const Frac &y) { return x -= y; }
	friend Frac operator * (Frac x, const Frac &y) { return x *= y; }
	friend Frac operator / (Frac x, const Frac &y) { return x /= y; }
	friend Frac operator - (const Frac &x) { return Frac(- x.x, x.y); }
	friend bool operator == (const Frac &x, const Frac &y) { return x.x * y.y == x.y * y.x; }
	friend bool operator != (const Frac &x, const Frac &y) { return x.x * y.y != x.y * y.x; }
	friend bool operator < (const Frac &x, const Frac &y) { return x.x * y.y < x.y * y.x; }
	friend bool operator > (const Frac &x, const Frac &y) { return x.x * y.y > x.y * y.x; }
	friend bool operator <= (const Frac &x, const Frac &y) { return x.x * y.y <= x.y * y.x; }
	friend bool operator >= (const Frac &x, const Frac &y) { return x.x * y.y >= x.y * y.x; }
	friend ostream& operator << (ostream& os, Frac x) { x.Reduce(); if (x.y == 1) return os << x.x; return os << x.x << '/' << x.y; }
	void Write() { Reduce(); if (y == 1) ::Write(x); else ::Write(x, '/', y); return; }
};
```