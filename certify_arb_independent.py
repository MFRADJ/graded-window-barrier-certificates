#!/usr/bin/env python3
"""Arb certificate for the three numerical regimes in Appendix A.

Requires python-flint==0.9.0. The domains are (0,4], [4,5], and
[5,infinity). All sign decisions use Arb balls and exit nonzero on
failure. The box bound is equation (A.3) of the accompanying preprint.
"""
from flint import arb, ctx
ctx.prec = 256
NH = 401                                      # half-step coefficients kept (exponent <= 200)

def mul(a, b, n):
    r = [0] * n
    for i, x in enumerate(a[:n]):
        if x:
            bi = b[:n - i]
            for j, y in enumerate(bi):
                r[i + j] += x * y
    return r

# ---- exact series in q^{1/2} ----------------------------------------------------------
th3 = [0] * NH; th4 = [0] * NH
for n in range(-30, 31):
    if n * n < NH:
        th3[n * n] += 1; th4[n * n] += 1 if n % 2 == 0 else -1
def pw4(a):
    s = mul(a, a, NH); return mul(s, s, NH)
T3, T4 = pw4(th3), pw4(th4)
T2 = [T3[i] - T4[i] for i in range(NH)]       # theta2^4 = theta3^4 - theta4^4 (checked exactly in series_exact.py)
M = (NH + 1) // 2
inv = [1] + [0] * (M - 1)
for n in range(1, M):
    for _ in range(4):
        for i in range(n, M):
            inv[i] += inv[i - n]
P4i = [0] * NH
for i in range(M):
    if 2 * i < NH: P4i[2 * i] = inv[i]
A2 = mul([(T3[i] + T4[i]) // 2 for i in range(NH)], P4i, NH)
B2 = mul([T2[i] // 2 for i in range(NH)], P4i, NH)
Ps = mul([T3[i] + T2[i] for i in range(NH)], P4i, NH)
a = [A2[2 * i] for i in range(M)]             # A = sum a_n q^n
c = [B2[2 * i + 1] for i in range(M - 1)]     # C = sum c_n q^n,  B = q^{1/2} C
if any(type(x) is not int for x in a + c + Ps):
    raise TypeError("Series coefficients must be exact integers")
if not (a[:4] == [1, 28, 134, 568] and c[:4] == [8, 64, 288, 1024]):
    raise AssertionError("Exact coefficient invariant failed")
if not (all(x >= 0 for x in a) and all(x >= 0 for x in c) and all(x >= 0 for x in Ps)):
    raise AssertionError("Exact coefficient invariant failed")

# ---- rigorous majorant constants at q0 = e^{-1} ------------------------------------------
q0 = arb(-1).exp()
r0 = q0.sqrt()
th3_up = 1 + 2 * r0 / (1 - r0)                       # sum_{n in Z} q0^{n^2/2} <= 1 + 2 sum_{k>=1} q0^{k/2}
th2_up = 2 * q0 ** (arb(1) / 8) / (1 - r0)           # (n+1/2)^2/2 >= n/2 + 1/8
invP4_up = (4 * q0 / (1 - q0) ** 2).exp()            # P(q0)^{-4}
GA = th3_up ** 4 * invP4_up                          # A(q0)   <= (theta3^4+theta4^4)/2 /P^4 <= theta3^4/P^4
GC = th3_up ** 4 * invP4_up / r0                     # C(q0)   <= q0^{-1/2} theta3^4 / P^4 (theta2^4 <= theta3^4)
GP = (th3_up ** 4 + th2_up ** 4) * invP4_up          # Psi(q0)

def tail_geo(r, N1):      # sum_{n>=N1} r^n
    return r ** N1 / (1 - r)
def tail_ngeo(r, N1):     # sum_{n>=N1} n r^n <= N1 r^N1/(1-r)^2  (valid for r<1, N1>=1)
    return N1 * r ** N1 / (1 - r) ** 2
def enclose(s, t):        # s + [0, t]
    return s + arb(0, 1) * t / 2 + t / 2

def A_delta_C(beta):
    """A(beta), delta=-dA/dbeta, C(e^-beta) as Arb balls; valid for beta >= 2."""
    beta = arb(beta); q = (-beta).exp(); r = q / q0
    N = len(a)
    A = sum(arb(a[n]) * q ** n for n in range(N))
    D = sum(arb(n * a[n]) * q ** n for n in range(1, N))
    C = sum(arb(c[n]) * q ** n for n in range(len(c)))
    A = enclose(A, GA * tail_geo(r, N)); D = enclose(D, GA * tail_ngeo(r, N))
    C = enclose(C, GC * tail_geo(r, len(c)))
    return A, D, C

def Sprime(y):
    """S'(y) = sum_{nu in Z/2, nu>=1/2} nu c_nu e^{-nu y}; valid for y >= 2."""
    y = arb(y); s = (-y / 2).exp(); r = s / r0          # half-step ratio
    S = sum(arb(i) / 2 * arb(Ps[i]) * s ** i for i in range(1, NH))
    return enclose(S, GP * tail_ngeo(r, NH) / 2)


def U(b1, b2, v1, v2):
    """Certified upper bound for the boxes, equation (A.3)."""
    (A1, D1, C1), (A2_, D2, C2) = v1, v2
    b1, b2 = arb(b1), arb(b2)
    return (A1.log() - b1 * D2 / A1
            + (1 + 2 * (-b1).exp() * C1 * (b2 * D1 / A2_).exp() / A2_).log())

def grid(lo, hi, n):
    return [arb(lo) + (arb(hi) - arb(lo)) * i / n for i in range(n + 1)]

if __name__ == "__main__":
    PI = arb.pi()
    print("Arb Appendix A certificate: (0,4] / [4,5] / [5,infty)")
    print("Precision: 256 bits; library: python-flint (Arb)")
    print("Majorants at q0=e^-1: GA = %s, GC = %s, GP = %s"
          % (GA.str(12), GC.str(12), GP.str(12)))
    y = PI ** 2
    R = y * (2 + y) * Sprime(y)
    lim = 4 * PI ** 2 / 3
    a_ok = bool(R < lim)
    print("(a) R(pi^2) = %s < 4pi^2/3 = %s : %s" % (R.str(12), lim.str(12), a_ok))
    P = grid(4, 5, 10)
    V = [A_delta_C(p) for p in P]
    ub = [U(P[i], P[i+1], V[i], V[i+1]) for i in range(10)]
    b_ok = all(u < 0 for u in ub)
    worst = max(u.upper() for u in ub)
    print("(b) [4,5], 10 boxes, bound (A.3): all < 0 : %s ; largest upper endpoint: %s"
          % (b_ok, worst.str(12)))
    A5, D5, C5 = V[-1]
    e5 = A5 - 1
    lhs = 28 * (4 - e5) / (1 + e5)
    rhs = 2 * C5 * (5 * D5).exp()
    c_ok = bool(lhs > rhs)
    print("(c) lhs = %s > rhs = %s; margin %s: %s" %
          (lhs.str(12), rhs.str(12), (lhs-rhs).str(12), c_ok))
    if not (a_ok and b_ok and c_ok):
        raise SystemExit("CERTIFICATE FAILED")
    print("==> APPENDIX A: CERTIFIED")
