#!/usr/bin/env python3
"""Contre-expertise S5.5 -- tests EXACTS (coefficients entiers) du theoreme 7.1
dans le regime mh >= 4000, jamais couvert par l'audit 6118 (m<=8, h<=150).
Variable x = q^{1/2}. Coefficients : entiers exacts (python-flint fmpz_poly).
Comparaisons finales : log2 flottant de grands entiers (test, non preuve),
lambda_h calcule exactement par bit_length. Sup sur TOUS les decalages t in Z/2
(Lemme 7.? : suffit), sans supposer t*=0.
Implementation independante : theta via produits de Jacobi, 1/P via pentagones."""
import math, sys
from flint import fmpz_poly

KAP = math.pi / math.log(2) * math.sqrt(8 / 3)
R2 = 1 / math.sqrt(2)

def trunc(p, L):
    return fmpz_poly([int(p[i]) for i in range(min(L, p.degree() + 1))])

def series(L):
    # 1/P(q) en x : P(x^2) = sum (-1)^k x^{k(3k-1)} (pentagonal), inverse par recurrence
    inv = [0] * L; inv[0] = 1
    pent = []
    k = 1
    while True:
        e1 = k * (3 * k - 1); e2 = k * (3 * k + 1)
        if e1 >= L: break
        s = 1 if k % 2 else -1
        pent.append((e1, s))
        if e2 < L: pent.append((e2, s))
        k += 1
    for n in range(2, L, 2):
        inv[n] = sum(s * inv[n - e] for e, s in pent if e <= n)
    IP = fmpz_poly(inv)
    IP4 = trunc(trunc(IP * IP, L) ** 2, L)
    # theta series (exposants x^{n^2})
    t3 = [0] * L; t4 = [0] * L; t2 = [0] * L
    n = 0
    while n * n < L:
        mult = 1 if n == 0 else 2
        t3[n * n] += mult; t4[n * n] += mult * (-1) ** n; n += 1
    T3 = fmpz_poly(t3); T4 = fmpz_poly(t4)
    T34 = trunc(trunc(T3 * T3, L) ** 2, L); T44 = trunc(trunc(T4 * T4, L) ** 2, L)
    num = T34 + T44
    A = trunc(num * IP4, L); A = fmpz_poly([int(A[i]) // 2 for i in range(L)])
    Bn = T34 - T44             # theta2^4 = theta3^4 - theta4^4 (Jacobi)
    B = trunc(Bn * IP4, L); B = fmpz_poly([int(B[i]) // 2 for i in range(L)])
    return A, B

def pw(p, e, L):
    r = fmpz_poly([1]); b = p
    while e:
        if e & 1: r = trunc(r * b, L)
        e >>= 1
        if e: b = trunc(b * b, L)
    return r

def lg(n): return math.log2(n) if n > 0 else 0.0

def run(m, enum, name, hs):
    hmax2 = max(hs)
    L = 4 * hmax2 + 2           # energie 4h -> indice 8h = 4*(2h)
    A, B = series(L)
    k = round(math.log2(sum(enum)))
    d = [0] * L
    for w, cw in enumerate(enum):
        if cw:
            ch = trunc(pw(A, m - w, L) * pw(B, w, L), L)
            for i in range(min(L, ch.degree() + 1)): d[i] += cw * int(ch[i])
    C = [0] * L; acc = 0
    for i in range(L): acc += d[i]; C[i] = acc
    def win(lo, hi):
        lo = max(lo, 0)
        if hi < lo: return 0
        hi = min(hi, L - 1)
        return C[hi] - (C[lo - 1] if lo > 0 else 0)
    out = []
    for h2 in hs:
        h = h2 / 2
        if d[h2] == 0:
            out.append((m, name, k, h, 'non admissible')); continue
        lam = d[h2].bit_length() - 1       # floor(log2) exact
        best = -1; tb = None
        for t2 in range(-2 * h2, 2 * h2 + 1):
            v = lg(max(1, win(t2, t2 + 2 * h2))) + lg(max(1, win(-t2, 2 * h2 - t2)))
            if v > best: best, tb = v, t2 / 2
        rho = best / (4 * lam)
        X = KAP * math.sqrt(m * h); e1 = 4 * m / 3 + 24 * (m * h) ** .25 + 1; e2 = e1 + 1
        low = R2 - 0.19 * math.sqrt(m / h) - 1.7 * (m * h) ** -.25
        up = R2 + 0.4 * math.sqrt(m / h) + 5 * (m * h) ** -.25
        upii = R2 * (X + k / math.sqrt(2)) / (X - e2) if X > e2 else float('inf')
        DS = lambda E2: lg(C[E2])
        # premisses (hd-level) et secteur-vacuum lower bound
        ok = (low <= rho <= up) and rho <= upii and lam >= X - e2 and DS(h2) <= X + k
        out.append((m, name, k, h, rho, rho - R2, low, up, upii, lam, round(X - e2, 1), tb, ok))
    return out

if __name__ == '__main__':
    cases = []
    comb = math.comb
    for m, hlist in ((1, [8000, 8001, 8002, 8003]), (2, [4000, 4001, 4002, 4003]),
                     (5, [1600, 1601]), (10, [800, 801]), (16, [500, 501])):
        encs = {'{0}': [1] + [0] * m,
                '<(1,0)>^m (k=m, impairs)': [comb(m, w) for w in range(m + 1)],
                'pair (k=m-1)': [comb(m, w) if w % 2 == 0 else 0 for w in range(m + 1)]}
        if m == 1: encs.pop('pair (k=m-1)')
        for nm, e in encs.items():
            cases += run(m, e, nm, hlist)
    bad = 0
    print("%3s %-26s %3s %8s %9s %10s %8s %8s %8s %6s %8s %6s %s" % ('m', 'code', 'k', 'h', 'rho', 'rho-1/r2', 'low(iii)', 'up(iii)', 'up(ii)', 'lam', 'X-e2', 't*', 'ok'))
    for r in cases:
        if r[4] == 'non admissible':
            print("%3d %-26s %3d %8.1f  non admissible" % r[:4]); continue
        m, nm, k, h, rho, dr, lo, up, u2, lam, xe, tb, ok = r
        bad += (not ok)
        print("%3d %-26s %3d %8.1f %9.6f %+10.6f %8.4f %8.4f %8.4f %6d %8.1f %6s %s" % (m, nm, k, h, rho, dr, lo, up, u2, lam, xe, tb, ok))
    print("cas testes:", sum(1 for r in cases if r[4] != 'non admissible'), " echecs:", bad)
    sys.exit(1 if bad else 0)
