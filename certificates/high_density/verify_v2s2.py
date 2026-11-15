#!/usr/bin/env python3
"""
Contre-expertise V2-S2 -- controles en arithmetique exacte.

Dependances : python-flint (entiers exacts), mpmath (controles flottants des
identites modulaires).  Execution :  python3 verify_v2s2.py > sortie.txt

Conventions (celles qui rendent coherentes la graduation demi-entiere et la
constante -log 2 du developpement modulaire) :
    q = e^{-beta},  x = q^{1/2},
    theta_3 = sum_n q^{n^2/2},  theta_4 = sum_n (-1)^n q^{n^2/2},
    theta_2 = sum_n q^{(n+1/2)^2/2},   P = prod_{n>=1} (1-q^n),
    A = (theta_3^4+theta_4^4)/(2 P^4),  B = theta_2^4/(2 P^4).
Toutes les series sont manipulees comme polynomes en x, tronques au degre X :
l'exposant n de x correspond au niveau n/2.
"""
import math
import sys
from flint import fmpz_poly
import mpmath as mp

LN2 = math.log(2.0)
a_const = 2 * math.pi ** 2 / 3            # a = 2 pi^2 / 3
KAPPA = (math.pi / LN2) * math.sqrt(8 / 3)  # kappa = (pi/ln 2) sqrt(8/3)
INV_SQRT2 = 1 / math.sqrt(2)


# --------------------------------------------------------------------------
# Series de base
# --------------------------------------------------------------------------
def base_series(X):
    """Renvoie (A, B) tronquees a x^X inclus (polynomes en x = q^{1/2})."""
    L = X + 1
    th3 = [0] * L
    n = 0
    while n * n < L:
        th3[n * n] += 1 if n == 0 else 2
        n += 1
    th3 = fmpz_poly(th3)
    th3_4 = th3.pow_trunc(4, L)
    c = th3_4.coeffs() + [0] * L
    # (theta3^4 + theta4^4)/2 = partie paire en x de theta3(x)^4
    N = fmpz_poly([int(c[i]) if i % 2 == 0 else 0 for i in range(L)])
    # theta2^4 / 2 = 8 x psi(x^2)^4,  psi(q) = sum_{n>=0} q^{n(n+1)/2}
    psi = [0] * L
    n = 0
    while n * (n + 1) < L:
        psi[n * (n + 1)] += 1
        n += 1
    M = (fmpz_poly(psi).pow_trunc(4, L) * 8).left_shift(1).truncate(L)
    # P(q)^{-1} = sum p(n) q^n, via la recurrence pentagonale ; en x : x^{2n}
    nq = X // 2 + 1
    p = [0] * nq
    p[0] = 1
    for i in range(1, nq):
        s = 0
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > i:
                break
            sign = 1 if k % 2 == 1 else -1
            s += sign * p[i - g1]
            if g2 <= i:
                s += sign * p[i - g2]
            k += 1
        p[i] = s
    Pinv = [0] * L
    for i in range(nq):
        if 2 * i < L:
            Pinv[2 * i] = p[i]
    Pinv4 = fmpz_poly(Pinv).pow_trunc(4, L)
    A = N.mul_low(Pinv4, L)
    B = M.mul_low(Pinv4, L)
    return A, B


def coeffs(poly, L):
    c = [int(v) for v in poly.coeffs()]
    return c + [0] * (L - len(c))


class Space:
    """Espace gradue V_S donne par son caractere (polynome en x)."""

    def __init__(self, name, m, k, Z, X):
        self.name, self.m, self.k, self.X = name, m, k, X
        self.d = coeffs(Z, X + 1)          # d[n] = dim (V_S)_{n/2}
        self.cum = [0] * (X + 1)
        s = 0
        for i, v in enumerate(self.d):
            s += v
            self.cum[i] = s

    def N(self, E):                        # dim (V_S)_{<=E}, E reel >= 0
        n = min(int(math.floor(2 * E + 1e-9)), self.X)
        return self.cum[n] if n >= 0 else 0

    def D(self, E):
        return math.log2(self.N(E))

    def dim_window(self, lo, hi):          # dim V[lo,hi], bornes en niveaux
        a = max(0, int(math.ceil(2 * lo - 1e-9)))
        b = int(math.floor(2 * hi + 1e-9))
        if b < a:
            return 0
        assert b <= self.X, "fenetre hors troncature"
        return self.cum[b] - (self.cum[a - 1] if a > 0 else 0)

    def W(self, lo, hi):
        return math.log2(max(1, self.dim_window(lo, hi)))

    def rho_win(self, h, lam=None):
        """Cout prescrit avec lambda maximal ; s parcourt (1/2)Z, |s|<=2h."""
        dh = self.d[int(round(2 * h))]
        if dh < 2:
            return None
        if lam is None:
            lam = dh.bit_length() - 1      # floor(log2 dim V_h)
        best, arg = -1.0, None
        for two_s in range(0, int(4 * h) + 1):   # symetrie s <-> -s
            s = two_s / 2
            v = self.W(s, s + 2 * h) + self.W(-s, 2 * h - s)
            if v > best:
                best, arg = v, s
        return best / (4 * lam), arg, lam


def build(name, m, k, terms, A, B, X):
    """terms = liste (multiplicite, i, j) : somme mult * A^i B^j."""
    L = X + 1
    Z = fmpz_poly([0])
    cacheA, cacheB = {}, {}
    for mult, i, j in terms:
        if i not in cacheA:
            cacheA[i] = A.pow_trunc(i, L)
        if j not in cacheB:
            cacheB[j] = B.pow_trunc(j, L)
        Z = Z + cacheA[i].mul_low(cacheB[j], L) * mult
    return Space(name, m, k, Z, X)


def binom_terms(m, base):
    """Developpe prod de blocs : base = liste de blocs, chacun une liste
    (mult, i, j).  Renvoie les termes du produit."""
    acc = {(0, 0): 1}
    for blk in base:
        new = {}
        for (i, j), c in acc.items():
            for mult, bi, bj in blk:
                key = (i + bi, j + bj)
                new[key] = new.get(key, 0) + c * mult
        acc = new
    return [(c, i, j) for (i, j), c in acc.items()]


TRIV = [(1, 1, 0)]                 # {0}           : A            (m=1,k=0)
GAM1 = [(1, 1, 0), (1, 0, 1)]      # {0,a}         : A+B          (m=1,k=1) auto-dual, impair
EPS2 = [(1, 2, 0), (3, 0, 2)]      # {00,aa,bb,cc} : A^2+3B^2     (m=2,k=2) auto-dual, pair


def code(name, blocks, A, B, X):
    m = sum(b[0][1] + b[0][2] for b in blocks)
    k = sum({1: 0, 2: 1}.get(len(b), 0) if b is not EPS2 else 2 for b in blocks)
    return build(name, m, k, binom_terms(m, blocks), A, B, X)


def macwilliams_selfdual(m, we):
    """Verifie W(x,y) = W((x+3y)/2,(x-y)/2) pour W = sum c_w x^(m-w) y^w
    (identite de MacWilliams d'un code kleinien auto-dual), en rationnels exacts."""
    from fractions import Fraction as Fr
    from math import comb
    out = [Fr(0)] * (m + 1)             # coefficients de x^(m-j) y^j
    for w, c in we.items():
        # ((x+3y)/2)^(m-w) * ((x-y)/2)^w
        p1 = [Fr(comb(m - w, i) * 3 ** i) for i in range(m - w + 1)]
        p2 = [Fr(comb(w, i) * (-1) ** i) for i in range(w + 1)]
        for i, u in enumerate(p1):
            for j, v in enumerate(p2):
                out[i + j] += c * u * v / 2 ** m
    return all(out[j] == we.get(j, 0) for j in range(m + 1))


# --------------------------------------------------------------------------
def section(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def check_basic(A, B):
    section("T1.a  Premiers coefficients (controle de normalisation)")
    L = 12
    cA, cB = coeffs(A, L), coeffs(B, L)
    print("A : coefficients de q^0..q^5      :", [cA[2 * i] for i in range(6)])
    print("B : coefficients de q^{1/2}..q^{9/2}:", [cB[2 * i + 1] for i in range(5)])
    print("A n'a que des puissances entieres de q :", all(cA[i] == 0 for i in range(1, L, 2)))
    E8 = coeffs(A.mul_low(A, L) + B.mul_low(B, L) * 3, L)
    print("A^2+3B^2 (doit etre le caractere de E8 : 1,248,4124,34752,213126) :",
          [E8[2 * i] for i in range(5)])


def check_modular():
    section("T1.b  Identites modulaires (mpmath, 40 chiffres)")
    mp.mp.dps = 40

    def AB(beta):
        q = mp.e ** (-beta)
        x = mp.sqrt(q)
        th3 = mp.nsum(lambda n: x ** (n * n), [-mp.inf, mp.inf])
        th4 = mp.nsum(lambda n: (-1) ** int(n) * x ** (n * n), [-mp.inf, mp.inf])
        th2 = mp.nsum(lambda n: x ** ((n + mp.mpf(1) / 2) ** 2), [-mp.inf, mp.inf])
        P = mp.qp(q)
        return (th3 ** 4 + th4 ** 4) / (2 * P ** 4), th2 ** 4 / (2 * P ** 4)

    a = 2 * mp.pi ** 2 / 3
    print(" beta | logA - [a/b - ln2 - b/6 + ln(A+3B)(q~)] | B/A - (A-B)/(A+3B)(q~) | reste exp. | a/b - logA")
    worst = 0
    for beta in [0.05, 0.1, 0.3, 0.7, 1, 2, 4, 2 * mp.pi, 8, 12, 20, 40]:
        A_, B_ = AB(beta)
        At, Bt = AB(4 * mp.pi ** 2 / beta)
        lhs = mp.log(A_)
        rhs = a / beta - mp.log(2) - beta / 6 + mp.log(At + 3 * Bt)
        r2 = B_ / A_ - (At - Bt) / (At + 3 * Bt)
        g = mp.log(At + 3 * Bt)
        worst = max(worst, abs(lhs - rhs), abs(r2))
        print(f"{float(beta):6.3f} | {mp.nstr(lhs - rhs, 3):>12} | {mp.nstr(r2, 3):>12} |"
              f" {mp.nstr(g, 5):>12} | {mp.nstr(a / beta - lhs, 8)}")
    print("ecart maximal des deux identites :", mp.nstr(worst, 3))
    # Lemme : log A(beta) <= a/beta pour tout beta>0 (grille fine)
    mn = min(float(a / b - mp.log(AB(b)[0])) for b in mp.linspace(0.02, 60, 1200))
    print("min sur une grille beta in [0.02,60] de a/beta - log A(beta) :", mn, "(>0 attendu)")
    # 0 < B/A < 1
    r = [float(AB(b)[1] / AB(b)[0]) for b in mp.linspace(0.02, 60, 300)]
    print("B/A sur la grille : min %.3e  max %.6f" % (min(r), max(r)))


def refined(sp, E):
    """Prediction heuristique : D ~ kappa sqrt(m(E-m/6)) - (m-k) - (1/4)log2(a m E) - log2(2 sqrt(pi))."""
    m, k = sp.m, sp.k
    Eeff = max(E - m / 6, 1e-9)
    return (KAPPA * math.sqrt(m * Eeff) - (m - k)
            - 0.25 * math.log2(a_const * m * E) - math.log2(2 * math.sqrt(math.pi)))


def check_V23(spaces):
    section("T2  Controle de (V2.3) : delta(E) = D_S(E) - kappa*sqrt(mE)")
    print("Colonnes : code | m | k=dim S | E | D_S(E) | delta | delta/(m+(mE)^(1/4)) | delta+(m-k) | D - prediction affinee")
    gmax_up, gmax_lo = -1e9, 1e9
    for sp in spaces:
        Emax = sp.X / 2
        Es = sorted(set([e for e in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 4000, 8000] if e <= Emax] + [Emax]))
        for E in Es:
            D = sp.D(E)
            delta = D - KAPPA * math.sqrt(sp.m * E)
            norm = delta / (sp.m + (sp.m * E) ** 0.25)
            print(f"{sp.name:>14} | {sp.m:3d} | {sp.k:3d} | {E:7.1f} | {D:10.2f} | {delta:9.3f} | {norm:7.3f} |"
                  f" {delta + sp.m - sp.k:8.3f} | {D - refined(sp, E):7.3f}")
        # extremum sur tous les E (entiers et demi-entiers)
        for n in range(0, sp.X + 1):
            E = n / 2
            delta = math.log2(sp.cum[n]) - KAPPA * math.sqrt(sp.m * E)
            gmax_up = max(gmax_up, delta - sp.k)        # borne sup. nette : D <= kappa sqrt(mE) + k
            gmax_lo = min(gmax_lo, delta / (sp.m + (sp.m * E) ** 0.25))
    print("\nSur tous les codes et tous les E testes :")
    print("  max de D_S(E) - kappa*sqrt(mE) - dim S          =", gmax_up, " (<=0 attendu : borne sup. nette)")
    print("  min de (D_S(E) - kappa*sqrt(mE))/(m+(mE)^(1/4)) =", gmax_lo, " (>= -C)")


def check_monotone(A, X):
    section("T3.a  Monotonie des coefficients de A^m (niveaux entiers)")
    L = X + 1
    for m in [1, 2, 3, 5, 8, 13]:
        c = coeffs(A.pow_trunc(m, L), L)
        ints = c[0::2]
        ok = all(ints[i] <= ints[i + 1] for i in range(len(ints) - 1))
        print(f"m={m:2d} : coefficients croissants jusqu'a q^{len(ints) - 1} : {ok}")


def check_windows(spaces, hs_of):
    section("T3.b  Cout des fenetres, encodage maximal lambda_h = floor(log2 dim V_h)")
    print("Colonnes : code | m | k | h | lambda_h | rho_win | rho_win-1/sqrt2 | argmax s | D(2h)/(2D(h)) - 1/sqrt2 | sqrt(m/h)")
    min_rho, min_ratio = 1e9, 1e9
    for sp in spaces:
        for h in hs_of(sp):
            if 4 * h > sp.X / 2:
                continue
            r = sp.rho_win(h)
            if r is None:
                print(f"{sp.name:>14} | {sp.m:3d} | {sp.k:3d} | {h:7.1f} | niveau non admissible (dim V_h < 2)")
                continue
            rho, arg, lam = r
            ratio = sp.D(2 * h) / (2 * sp.D(h))
            min_rho = min(min_rho, rho)
            min_ratio = min(min_ratio, ratio)
            print(f"{sp.name:>14} | {sp.m:3d} | {sp.k:3d} | {h:7.1f} | {lam:6d} | {rho:.5f} | {rho - INV_SQRT2:+.5f} |"
                  f" {arg:6.1f} | {ratio - INV_SQRT2:+.5f} | {math.sqrt(sp.m / h):.4f}")
    print("\nmin de rho_win observe             :", min_rho, "  (1/sqrt2 = %.6f)" % INV_SQRT2)
    print("min de D(2h)/(2D(h)) observe       :", min_ratio)


def exhaustive_from_above(spaces):
    section("T4  Recherche systematique : rho_win < 1/sqrt2 ou D(2h) < sqrt2*D(h) ?")
    bad_rho, bad_ratio, tested = [], [], 0
    for sp in spaces:
        hmax = sp.X / 8
        n = 2
        while n / 2 <= hmax:
            h = n / 2
            n += 1
            if sp.d[int(2 * h)] < 2:
                continue
            tested += 1
            if sp.D(2 * h) < math.sqrt(2) * sp.D(h):
                bad_ratio.append((sp.name, h))
            # rho_win >= fenetre centrale / (4 lambda) : test complet seulement si necessaire
            lam = sp.d[int(2 * h)].bit_length() - 1
            if 2 * sp.D(2 * h) / (4 * lam) < INV_SQRT2:
                r = sp.rho_win(h)
                if r[0] < INV_SQRT2:
                    bad_rho.append((sp.name, h, r[0]))
    print("couples (code, h) admissibles testes :", tested)
    print("violations de D(2h) >= sqrt2*D(h)    :", len(bad_ratio), bad_ratio[:10])
    print("violations de rho_win >= 1/sqrt2     :", len(bad_rho), bad_rho[:10])


def main():
    X1 = 2 * 8000                       # m = 1 : niveaux jusqu'a 8000
    A1, B1 = base_series(X1)
    check_basic(A1, B1)
    check_modular()

    sp = []
    sp.append(code("trivial", [TRIV], A1, B1, X1))
    sp.append(code("gamma1", [GAM1], A1, B1, X1))
    X = 2 * 4000
    A, B = A1.truncate(X + 1), B1.truncate(X + 1)
    sp.append(code("trivial", [TRIV] * 2, A, B, X))
    sp.append(code("eps2(E8)", [EPS2], A, B, X))
    sp.append(code("gamma1^2", [GAM1] * 2, A, B, X))
    sp.append(code("trivial", [TRIV] * 4, A, B, X))
    sp.append(code("eps2^2", [EPS2] * 2, A, B, X))
    sp.append(code("eps2+triv^2", [EPS2] + [TRIV] * 2, A, B, X))
    sp.append(code("gamma1^3+triv", [GAM1] * 3 + [TRIV], A, B, X))
    X = 2 * 2000
    A, B = A1.truncate(X + 1), B1.truncate(X + 1)
    for m in (8, 12):
        sp.append(code("trivial", [TRIV] * m, A, B, X))
        sp.append(code("eps2^%d" % (m // 2), [EPS2] * (m // 2), A, B, X))
        sp.append(code("gamma1^%d" % m, [GAM1] * m, A, B, X))
        sp.append(code("mixte k=m/2", [EPS2] * (m // 4) + [TRIV] * (m // 2), A, B, X))
    # Codes auto-duaux pairs indecomposables, donnes par leur enumerateur de poids
    # (invariance de MacWilliams verifiee exactement ci-dessous).
    hexa = {0: 1, 4: 45, 6: 18}                       # hexacode, m = 6
    dodeca = {0: 1, 6: 396, 8: 1485, 10: 1980, 12: 234}  # dodecacode, m = 12
    for nm, mm, we in (("hexacode", 6, hexa), ("dodecacode", 12, dodeca)):
        assert macwilliams_selfdual(mm, we), nm
        sp.append(build(nm, mm, mm, [(c, mm - w, w) for w, c in we.items()], A, B, X))
    X = 2 * 3600                        # m = 30, h = m^2 = 900, fenetres jusqu'a 4h
    A, B = A1.truncate(X + 1), B1.truncate(X + 1)
    for m in (20, 30):
        sp.append(code("trivial", [TRIV] * m, A, B, X))
        sp.append(code("eps2^%d" % (m // 2), [EPS2] * (m // 2), A, B, X))
        sp.append(code("gamma1^%d" % m, [GAM1] * m, A, B, X))

    check_V23(sp)
    check_monotone(A1, 2 * 1500)

    def hs_of(s):
        m = s.m
        cand = {m, 2 * m, 5 * m, 10 * m, 50 * m, m * m, round(m * math.log(max(m, 2))),
                10, 100, 500, 1000, 2000, 100.5, 500.5}
        return sorted(h for h in cand if h >= 1)

    check_windows(sp, hs_of)
    exhaustive_from_above(sp)


if __name__ == "__main__":
    main()
