```cpp
vector<int> minp, primes, lowp;

void Sieve(int n) {
	minp.assign(n, 0);
	lowp.assign(n, 0);
	primes.clear();
	
	for (int i = 2; i < n; i ++) {
		if (minp[i] == 0) {
			minp[i] = i, lowp[i] = i;
			primes.emplace_back(i);
		}
		
		for (auto p : primes) {
			if (i * p >= n) {
				break;
			}
			minp[i * p] = p;
			if (p == minp[i]) {
				lowp[i * p] = lowp[i] * p;
				break;
			}
			lowp[i * p] = p;
		}
	}
}
bool is_prime(int n) {
	return minp[n] == n;
}
```