对于一个矩阵 $\mathbf A$，可以使用 `charPoly(A)` 在 $O(n^3)$ 的时间复杂度内求出 $\mathbf A$ 的特征多项式 $\det(\lambda \mathbf I - \mathbf A)$。

对于 $m + 1$ 个矩阵 $\mathbf A_0, \cdots, \mathbf A_m$，可以使用 `detPoly(A)` 在 $O(n^3m^3)$ 的时间复杂度内求出 $\det(\mathbf A_0 + \mathbf A_1 x + \mathbf A_2 x^2 + \cdots + \mathbf A_m x^m)$。

```cpp
std::vector<Z> charPoly(std::vector<std::vector<Z>> A) {
	const size_t n = A.size();
	for (size_t i = 0; i < n; i ++) {
		for (size_t j = 0; j < n; j ++) {
			A[i][j] = -A[i][j];
		}
	}
	for (size_t i = 0; i + 2 < n; i ++) {
		size_t pivot = i + 1;
		for (; pivot < n && !A[pivot][i]; ++ pivot);
		if (pivot == n) {
			continue;
		}
		if (pivot > i + 1) {
			for (size_t j = i; j < n; j ++) {
				std::swap(A[i + 1][j], A[pivot][j]);
			}
			for (size_t j = 0; j < n; j ++) {
				std::swap(A[j][i + 1], A[j][pivot]);
			}
		}
		const Z inv = A[i + 1][i].inv();
		for (size_t j = i + 2; j < n; j ++) {
			if (A[j][i]) {
				const Z t = A[j][i] * inv;
				for (size_t k = i; k < n; k ++) {
					A[j][k] = fsm(A[j][k], t, A[i + 1][k]);
				}
				for (size_t k = 0; k < n; k ++) {
					A[k][i + 1] = fam(A[k][i + 1], t, A[k][j]);
				}
			}
		}
	}
	std::vector<std::vector<Z>> dp(n + 1);
	dp[0] = {Z(1)};
	for (size_t i = 0; i < n; i ++) {
		dp[i + 1].assign(i + 2, Z(0));
		for (size_t k = 0; k <= i; k ++) {
			dp[i + 1][k + 1] = dp[i][k];
		}
		for (size_t k = 0; k <= i; k ++) {
			dp[i + 1][k] = fam(dp[i + 1][k], A[i][i], dp[i][k]);
		}
		Z prod = 1;
		for (size_t j = i; j; j --) {
			prod *= -A[j][j - 1];
			const Z t = prod * A[j - 1][i];
			for (size_t k = 0; k < j; k ++) {
				dp[i + 1][k] = fam(dp[i + 1][k], t, dp[j - 1][k]);
			}
		}
	}
	return dp[n];
}
std::vector<Z> detPoly(std::vector<std::vector<std::vector<Z>>> A) {
	assert(A.size());
	const int m = A.size() - 1, n = A[0].size();
	Z prod = Z(1); int offset = 0;
	for (int h = 0; h < n; h ++) {
		for (; ; ) {
			if (A[m][h][h]) {
				break;
			}
			for (int j = h + 1; j < n; j ++) {
				if (A[m][h][j]) {
					prod = -prod;
					for (int d = 0; d <= m; d ++) {
						for (int i = 0; i < n; i ++) {
							std::swap(A[d][i][h], A[d][i][j]);
							break;
						}
					}
				}
			}
			if (A[m][h][h]) {
				break;
			}
			if (++ offset > n * m) {
				return std::vector<Z>(n * m + 1, Z(0));
			}
			for (int d = m; -- d >= 0; ) {
				for (int j = 0; j < n; j ++) {
					A[d + 1][h][j] = A[d][h][j];
				}
			}
			for (int j = 0; j < n; j ++) {
				A[0][h][j] = Z(0);
			}
			for (int i = 0; i < h; i ++) {
				const Z t = A[m][h][i];
				for (int d = 0; d <= m; d ++) {
					for (int j = 0; j < n; j ++) {
						A[d][h][j] -= t * A[d][i][j];
					}
				}
			}
		}
		prod *= A[m][h][h];
		const Z inv = A[m][h][h].inv();
		for (int d = 0; d <= m; d ++) {
			for (int j = 0; j < n; j ++) {
				A[d][h][j] *= inv;
			}
		}
		for (int i = 0; i < n; i ++) {
			if (h != i) {
				const Z t = A[m][i][h];
				for (int d = 0; d <= m; d ++) {
					for (int j = 0; j < n; j ++) {
						A[d][i][j] -= t * A[d][h][j];
					}
				}
			}
		}
	}
	std::vector<std::vector<Z>> B(n * m, std::vector<Z>(n * m, Z(0)));
	for (int i = 0; i < (m - 1) * n; i ++) {
		B[i][n + i] = Z(1);
	}
	for (int d = 0; d < m; d ++) {
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < n; j ++) {
				B[(m - 1) * n + i][d * n + j] = -A[d][i][j];
			}
		}
	}
	const std::vector<Z> f = charPoly(B);
	std::vector<Z> g(n * m + 1, Z(0));
	for (int i = 0; i + offset <= n * m; i ++) {
		g[i] = prod * f[i + offset];
	}
	return g;
}
```