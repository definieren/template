```cpp
Z Pfaffian(vector<vector<Z>> A) {
	const int n = A.size();
	assert((n & 1) == 0);
	
	bool sgn = false;
	Z Pf = 1;
	for (int i = 0; i < n; i += 2) {
		int p = -1;
		for (int j = i + 1; j < n; j ++) {
			if (A[i][j]) {
				p = j;
				break;
			}
		}
		if (p == -1) {
			return 0;
		}
		if (p != i + 1) {
			sgn ^= 1;
			for (int j = 0; j < n; j ++) {
				swap(A[j][i + 1], A[j][p]);
			}
			for (int j = i + 1; j < n; j ++) {
				swap(A[i + 1][j], A[p][j]);
			}
		}
		const Z inv = A[i][i + 1].inv();
		for (int j = i + 2; j < n; j ++) {
			if (A[i][j]) {
				const Z coef = -A[i][j] * inv;
				for (int k = i; k < n; k ++) {
					A[k][j] += A[k][i + 1] * coef;
				}
				for (int k = 0; k < n; k ++) {
					A[j][k] += A[i + 1][k] * coef;
				}
			}
		}
		Pf *= A[i][i + 1];
	}
	return sgn ? -Pf : Pf;
}
```