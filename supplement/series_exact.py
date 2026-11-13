"""Exact integer checks of the q-series identities used in Appendix A.

Series are stored in powers of q^{1/2}: index i represents exponent i/2.
The checks cover every exponent strictly below 200."""
N = 400                     # number of half-units kept (exponent < 200)

def check(description, passed):
    print(description, passed)
    if not passed:
        raise SystemExit("EXACT SERIES CHECK FAILED")
def mul(a, b):
    r = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                r[i+j] += x*b[j]
    return r
def power(a, k):
    r = [1] + [0]*(N-1)
    for _ in range(k): r = mul(r, a)
    return r
# theta3 = sum q^{n^2/2} -> index n^2 ; theta4 with (-1)^n
th3 = [0]*N; th4 = [0]*N
n = 0
while n*n < N:
    for s in ([n] if n == 0 else [n, -n]):
        th3[n*n] += 1; th4[n*n] += (-1)**n
    n += 1
# theta2 = sum q^{(n+1/2)^2/2}: exponent (2n+1)^2/8 is not a multiple of 1/2 -> compute theta2^4 directly by 4-fold sum
th2_4 = [0]*N
odd = [k for k in range(-60, 61) if k % 2]   # k = 2n+1
from itertools import product
sq = {k: k*k for k in odd}
for a, b, c, d in product(odd, repeat=4):
    tot = sq[a]+sq[b]+sq[c]+sq[d]          # = 8 * exponent ; exponent*2 = tot/4
    if tot % 4: raise SystemExit("parity")
    i = tot // 4
    if i < N: th2_4[i] += 1
th3_4 = power(th3, 4); th4_4 = power(th4, 4)
jac = [th3_4[i] - th4_4[i] - th2_4[i] for i in range(N)]
check("Jacobi quartic theta3^4 = theta2^4 + theta4^4 for exponents < 200:", all(x == 0 for x in jac))
# 1/P^4 in powers of q (integer exponents) -> place at even indices
M = N//2
inv = [1] + [0]*(M-1)
for n_ in range(1, M):          # multiply by 1/(1-q^n)^4
    for _ in range(4):
        for i in range(n_, M): inv[i] += inv[i-n_]
P4inv = [0]*N
for i in range(M): P4inv[2*i] = inv[i]
Anum = [(th3_4[i] + th4_4[i])//2 for i in range(N)]
Bnum = [(th3_4[i] - th4_4[i])//2 for i in range(N)]
if not (all((th3_4[i] + th4_4[i]) % 2 == 0 for i in range(N))):
    raise AssertionError("Exact coefficient invariant failed")
A = mul(Anum, P4inv); B = mul(Bnum, P4inv)
check("A only integer powers of q :", all(A[i] == 0 for i in range(1, N, 2)))
check("B only half-odd powers     :", all(B[i] == 0 for i in range(0, N, 2)))
Aq = [A[2*i] for i in range(M)]; Cq = [B[2*i+1] for i in range(M)]
print("A first coeffs:", Aq[:6]); print("C first coeffs:", Cq[:6])
check("A and C coefficients are exact integers:", all(type(x) is int for x in Aq + Cq))
check("A first coefficients match:", Aq[:6] == [1, 28, 134, 568, 1809, 5316])
check("C first coefficients match:", Cq[:6] == [8, 64, 288, 1024, 3152, 8704])
check("A-1 >= 0 and C >= 0 up to q^%d:" % (M-1), all(x >= 0 for x in Aq[1:]) and all(x >= 0 for x in Cq))
Psi = mul([th3_4[i] + th2_4[i] for i in range(N)], P4inv)
print("Psi first coeffs (half-steps):", Psi[:8])
check("Psi coefficients are nonnegative integers:", all(type(x) is int and x >= 0 for x in Psi))
check("A+3B == Psi :", all(A[i] + 3*B[i] == Psi[i] for i in range(N)))
# Jacobi triple product for theta3, theta4 (q^{1/2} variable): prod (1-q^n)(1 +/- q^{n-1/2})^2
def tp(sign):
    r = [1] + [0]*(N-1)
    for n_ in range(1, N):
        f = [0]*N; f[0] = 1
        if 2*n_ < N: f[2*n_] = -1
        r = mul(r, f)
        if 2*n_-1 < N:
            g = [0]*N; g[0] = 1; g[2*n_-1] = sign
            r = mul(r, g); r = mul(r, g)
    return r
check("Triple product theta3:", tp(+1) == th3)
check("Triple product theta4:", tp(-1) == th4)

print("==> EXACT SERIES CHECKS: PASS")
