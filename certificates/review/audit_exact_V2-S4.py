#!/usr/bin/env python3
"""V2-S4 -- audit exact (entiers) du theoreme V2-S3.
Variable x = q^{1/2} ; l'indice n represente l'energie n/2.
A = (th3^4+th4^4)/(2P^4), B = th2^4/(2P^4) (conventions S2-bis), recalcules ici
directement par produits de series, sans reutiliser le code du dossier.
dim(V_S)_{n/2} = sum_w c_w(S) [x^n] A^{m-w} B^w, c = enumerateur de poids de S.
"""
import itertools, math, random, sys
from flint import fmpz_poly

KAPPA = (math.pi / math.log(2)) * math.sqrt(8 / 3)
ISQ2 = 1 / math.sqrt(2)


def series(L):
    """A, B en x (x = q^{1/2}) tronques a l'ordre L (indices 0..L-1)."""
    th3 = [0] * L; th4 = [0] * L; psi = [0] * L
    n = 0
    while n * n < L:
        for s in ((n, -n) if n else (0,)):
            th3[n * n] += 1; th4[n * n] += (-1) ** n
        n += 1
    n = 0
    while n * n + n < L:          # th2 = x^{1/4}... : th2^4 = x * (sum_{n in Z} x^{n^2+n})^4
        psi[n * n + n] += 2       # n et -n-1 donnent le meme exposant
        n += 1
    T3, T4, Ps = fmpz_poly(th3), fmpz_poly(th4), fmpz_poly(psi)
    # 1/P^4 avec P = prod (1 - x^{2n}) : on calcule les partitions 4-colorees en q
    inv = [0] * L; inv[0] = 1
    for _ in range(4):
        for k in range(2, L, 2):
            for i in range(k, L):
                inv[i] += inv[i - k]
    IP = fmpz_poly(inv)
    num = (T3 ** 4 + T4 ** 4)
    A = (num * IP)
    A = fmpz_poly([int(A[i]) // 2 for i in range(L)])
    th2_4 = fmpz_poly([0] + [int((Ps ** 4)[i]) for i in range(L - 1)])  # x * psi^4
    B = th2_4 * IP
    B = fmpz_poly([int(B[i]) // 2 for i in range(L)])
    return A, B


def trunc(p, L):
    return fmpz_poly([int(p[i]) for i in range(min(L, p.degree() + 1))])


def mul(p, q, L):
    return trunc(p * q, L)


def powtrunc(p, e, L):
    r = fmpz_poly([1]); b = p
    while e:
        if e & 1: r = mul(r, b, L)
        e >>= 1
        if e: b = mul(b, b, L)
    return r


def sector_chars(m, L, A, B):
    """liste [A^{m-w} B^w pour w=0..m] tronquee a L."""
    Ap = [fmpz_poly([1])]
    for _ in range(m): Ap.append(mul(Ap[-1], A, L))
    Bp = [fmpz_poly([1])]
    for _ in range(m): Bp.append(mul(Bp[-1], B, L))
    return [mul(Ap[m - w], Bp[w], L) for w in range(m + 1)]


def coeffs(p, L):
    return [int(p[i]) if i <= p.degree() else 0 for i in range(L)]


def lg(x):
    return math.log2(x) if x > 0 else float('-inf')


# ---------- enumerateurs de poids des sous-groupes de K^m ----------
def wt(v, m):  # v entier sur 2m bits, coordonnee i = bits (2i, 2i+1)
    return sum(1 for i in range(m) if (v >> (2 * i)) & 3)


def span(gens):
    S = {0}
    for g in gens:
        S |= {s ^ g for s in S}
    return S


def all_subgroups(m, kmax):
    """tous les sous-groupes de F_2^{2m} de dim <= kmax (via formes echelonnees)."""
    seen = set(); out = []
    vecs = range(1, 1 << (2 * m))
    out.append(frozenset({0}))
    for k in range(1, kmax + 1):
        for gens in itertools.combinations(vecs, k):
            S = frozenset(span(gens))
            if len(S) != 1 << k or S in seen: continue
            seen.add(S); out.append(S)
    return out


def enumerator(S, m):
    c = [0] * (m + 1)
    for s in S: c[wt(s, m)] += 1
    return tuple(c)


def symp(u, v, m):  # forme alternee b(x,y)=x1y2+x2y1 par coordonnee
    t = 0
    for i in range(m):
        a1, a2 = (u >> 2 * i) & 1, (u >> 2 * i + 1) & 1
        b1, b2 = (v >> 2 * i) & 1, (v >> 2 * i + 1) & 1
        t ^= (a1 & b2) ^ (a2 & b1)
    return t


# ---------- audit d'un couple (m, enumerateur) ----------
def audit(m, c, chars, Lidx, hs, log):
    k = int(round(math.log2(sum(c))))
    d = [0] * Lidx
    for w, cw in enumerate(c):
        if cw:
            cf = coeffs(chars[w], Lidx)
            for i in range(Lidx): d[i] += cw * cf[i]
    C = list(itertools.accumulate(d))
    secC = {w: list(itertools.accumulate(coeffs(chars[w], Lidx))) for w in range(m + 1) if c[w]}
    worst = {}
    def upd(key, val, info):
        if key not in worst or val < worst[key][0]: worst[key] = (val, info)
    viol = []
    for h2 in hs:                       # h2 = 2h (indice de niveau)
        h = h2 / 2
        if 4 * h2 >= Lidx: break        # fenetres jusqu'a l'energie 4h
        if d[h2] == 0:
            continue                    # niveau non admissible
        X = KAPPA * math.sqrt(m * h)
        e1 = 4 * m / 3 + 24 * (m * h) ** 0.25 + 1; e2 = e1 + 1
        lam = int(math.floor(lg(d[h2])))
        # exactitude du floor sur les grands entiers
        while (1 << (lam + 1)) <= d[h2]: lam += 1
        while (1 << lam) > d[h2]: lam -= 1
        DS = lambda E2: lg(C[E2])       # E2 = 2E
        # entrees (1)-(3)
        for E2 in (h2, 2 * h2):
            E = E2 / 2
            up = KAPPA * math.sqrt(m * E) + k - DS(E2); upd('(1) D_S <= kappa sqrt(mE)+k', up, (m, c, E))
            if up < -1e-9: viol.append(('(1)', m, c, E, up))
            for w in secC:               # borne par secteur (Prop. 3.2 S2-bis)
                lo = lg(secC[w][E2]) - (KAPPA * math.sqrt(m * E) - 4 * m / 3 - 20 * (m * E) ** 0.25 - 1)
                upd('(2) log2 N_s >= ...', lo, (m, c, E, w))
                if lo < -1e-9: viol.append(('(2)', m, c, E, w, lo))
        g = lam - (X - e2); upd('(3) lambda_h >= X - e2', g, (m, c, h))
        if g < -1e-9: viol.append(('(3a)', m, c, h, g))
        g = (X + k) - DS(h2); upd('(3) D_S(h) <= X+k', g, (m, c, h))
        # rho_win exact (sup sur t dans Z/2, Lemme 4.3 S2-bis)
        def win(lo2, hi2):               # dim V[lo2/2, hi2/2]
            lo2 = max(lo2, 0)
            if hi2 < lo2: return 0
            return C[hi2] - (C[lo2 - 1] if lo2 > 0 else 0)
        best = -1; tbest = None
        for t2 in range(-2 * h2, 2 * h2 + 1):
            v = lg(max(1, win(t2, t2 + 2 * h2))) + lg(max(1, win(-t2, 2 * h2 - t2)))
            if v > best: best, tbest = v, t2 / 2
        rho = best / (4 * lam)
        lowI = DS(2 * h2) / (2 * DS(h2))
        lowIb = ISQ2 * (X - e1 / math.sqrt(2)) / (X + k)
        lowIc = ISQ2 - 0.19 * math.sqrt(m / h) - 1.7 * (m * h) ** -0.25
        chain = [rho, lowI, lowIb, lowIc]
        for j in range(3):
            if chain[j] < chain[j + 1] - 1e-12: viol.append(('(i) chaine', j, m, c, h, chain))
        if X > e2:
            upII = ISQ2 * (X + k / math.sqrt(2)) / (X - e2)
            if rho > upII + 1e-12: viol.append(('(ii)', m, c, h, rho, upII))
        if m * h >= 4000:
            upIII = ISQ2 + 0.4 * math.sqrt(m / h) + 5 * (m * h) ** -0.25
            if not (lowIc - 1e-12 <= rho <= upIII + 1e-12): viol.append(('(iii)', m, c, h, rho))
        upd('rho - 1/sqrt2', rho - ISQ2, (m, c, h, 'demi' if h2 % 2 else 'entier', 't*=%s' % tbest))
        log.append((m, c, k, h, rho, tbest, lam, X))
    return worst, viol


def main():
    A, B = series(64)
    # controle des premiers coefficients (S2-bis) : A = 1+28q+134q^2, B = 8q^{1/2}(1+8q+36q^2)
    a_ = coeffs(A, 8); b_ = coeffs(B, 8)
    assert a_[:6] == [1, 0, 28, 0, 134, 0], a_
    assert b_[:7] == [0, 8, 0, 64, 0, 288, 0], b_
    # A - B = th4^4/P^4 > 0 et monotonie par classe : controles rapides
    print('A :', a_, ' B :', b_)
    allworst = {}; allviol = []; LOG = []
    plan = []
    # (a) exhaustif m=1,2,3 (k<=m), h de 15m a Hmax
    for m, Hmax in ((1, 150), (2, 100), (3, 80)):
        subs = all_subgroups(m, m)
        enums = sorted({enumerator(S, m) for S in subs})
        iso = sum(1 for S in subs if all(symp(u, v, m) == 0 for u in S for v in S))
        print('m=%d : %d sous-groupes (k<=m), %d isotropes, %d enumerateurs distincts' % (m, len(subs), iso, len(enums)))
        plan.append((m, enums, Hmax))
    # (b) m = 4..8, enumerateurs aleatoires + extremes
    rnd = random.Random(20261114)
    for m, Hmax in ((4, 100), (5, 110), (6, 120), (8, 150)):
        enums = set()
        enums.add(tuple([1] + [0] * m))                                   # {0}
        enums.add(tuple(math.comb(m, w) for w in range(m + 1)))           # <(1,0)_i> : k=m, mots impairs
        enums.add(tuple((math.comb(m, w) if w % 2 == 0 else 0) for w in range(m + 1)))  # pairs, k=m-1
        for _ in range(12):
            k = rnd.randint(1, m)
            gens = [rnd.randrange(1, 1 << (2 * m)) for _ in range(k)]
            enums.add(enumerator(span(gens), m))
        plan.append((m, sorted(enums), Hmax))
    for m, enums, Hmax in plan:
        Lidx = 8 * Hmax + 2
        A, B = series(Lidx)
        chars = sector_chars(m, Lidx, A, B)
        hs = range(2 * 15 * m, 2 * Hmax + 1)
        for c in enums:
            w_, v_ = audit(m, c, chars, Lidx, hs, LOG)
            allviol += v_
            for key, val in w_.items():
                if key not in allworst or val[0] < allworst[key][0]: allworst[key] = val
        print('  m=%d termine (%d enumerateurs, h jusqu a %d)' % (m, len(enums), Hmax))
    print('\n== Marges minimales observees (>= 0 attendu)')
    for key, (val, info) in allworst.items():
        print('  %-32s min = %+.4f  en %s' % (key, val, info))
    print('\n== Violations : %d' % len(allviol))
    for v in allviol[:20]: print('  ', v)
    print('couples (code, niveau) testes :', len(LOG))
    below = [r for r in LOG if r[4] < ISQ2]
    print('cas rho < 1/sqrt2 :', len(below))
    tnon0 = [r for r in LOG if r[5] != 0]
    print('cas ou le sup n est PAS atteint en t=0 :', len(tnon0), tnon0[:5])
    with open('audit_log.tsv', 'w') as f:
        f.write('m\tenum\tk\th\trho\tt*\tlambda_h\tX\n')
        for r in LOG: f.write('\t'.join(map(str, r)) + '\n')
    return 1 if allviol else 0


if __name__ == '__main__':
    sys.exit(main())
