#!/usr/bin/env python3
"""Contre-expertise V2-S2-bis -- certificat Arb (python-flint) de toutes les constantes
decimales du manuscrit, y compris celles que certify_bis1_bis2.py ne controle pas.
Chaque test est une comparaison de boules Arb (arithmetique d'intervalles rigoureuse) :
un test PASS signifie que l'inegalite est prouvee, pas seulement observee.
Deuxieme partie : controle numerique (mpmath, 60 chiffres) des quantites modulaires
du lemme 2.1 -- informatif seulement, la preuve du lemme etant analytique.
"""
import sys
from flint import arb, ctx
import mpmath as mp

ctx.prec = 200
res = []


def chk(name, ok):
    ok = bool(ok)
    res.append(ok)
    print(('PASS ' if ok else 'FAIL ') + name)


def lt(x, y):  # x < y prouve
    return (arb(y) - arb(x)) > 0


pi = arb.pi()
ln2 = arb(2).log()
eps = arb(1) / 10000
a = 2 * pi ** 2 / 3
c = 2 * pi ** 2
kappa = 2 * a.sqrt() / ln2
sq2 = arb(2).sqrt()

print('== Lemme 2.1 : constantes de y, y\', y\'\'')
chk('19.7 < c=2pi^2 < 20', lt(arb('19.7'), c) and lt(c, 20))
chk('e^{-c} < 3e-9', lt((-c).exp(), arb(3) / 10 ** 9))
chk("y' <= c e^{-c} < 6e-8", lt(c * (-c).exp(), arb(6) / 10 ** 8))
chk("|y''| <= c^2 e^{-c} < 1.32e-6  (non controle par le script du dossier)", lt(c ** 2 * (-c).exp(), arb('1.32e-6')))
chk('c/4 > 1 (monotonie de beta^{-j} e^{-c/beta} sur (0,1], j<=4)', lt(4, c))
chk('c > 2 (|c^2/b^4 - 2c/b^3| <= c^2/b^4 sur (0,1])', lt(2, c))
R = arb(1) / 16
v = R * R
P4 = 1 / (1 - 4 * v / (1 - v) ** 2)
Abd = ((1 + R) / (1 - R)) ** 4 * P4
Bbd = 8 * R / (1 - v) ** 4 * P4
chk('|A(z^2)| <= 1.68 sur |z|=1/16', lt(Abd, arb('1.68')))
chk('|B(z^2)| <= 0.52 sur |z|=1/16', lt(Bbd, arb('0.52')))
chk('|F|,|G| < 4 sur |z|=1/16', lt(Abd + 3 * Bbd, 4))
print('== Lemme 2.1 : borne globale A < e^{a/beta}')
Rs = (-pi).exp()
qs = (-2 * pi).exp()
chk('e^{-pi} < 1/23', lt(Rs, arb(1) / 23))
chk('e^{-2pi} < 1/500', lt(qs, arb(1) / 500))
Aqs = (((1 + Rs) / (1 - Rs)) ** 4 + ((1 - Rs) / (1 + Rs)) ** 4) / 2 / (1 - 4 * qs / (1 - qs) ** 2)
Bqs = 8 * Rs / (1 - qs) ** 4 / (1 - 4 * qs / (1 - qs) ** 2)
chk('majorant de A(e^{-2pi}) < 1.08 (valeur %s)' % Aqs.str(6), lt(Aqs, arb('1.08')))
chk('majorant de B(e^{-2pi}) < 0.36 (valeur %s)' % Bqs.str(6), lt(Bqs, arb('0.36')))
chk('F(e^{-pi}) < 2.2 < 2 e^{1/6}  (2e^{1/6} non controle par le script du dossier)',
    lt(Aqs + 3 * Bqs, arb('2.2')) and lt(arb('2.2'), 2 * (arb(1) / 6).exp()))
chk('0.08 * 2pi < a', lt(arb('0.08') * 2 * pi, a))
print('== Proposition 3.2')
chk('phi(1)/m <= a + 1/6 + eps + 2 sqrt(2a+2eps) < 14.01 (valeur %s)' % (a + arb(1) / 6 + eps + 2 * (2 * a + 2 * eps).sqrt()).str(8),
    lt(a + arb(1) / 6 + eps + 2 * (2 * a + 2 * eps).sqrt(), arb('14.01')))
chk('1 + 2eps/ln2 < 1.001', lt(1 + 2 * eps / ln2, arb('1.001')))
chk('(1/3+2eps) sqrt(a)/ln2 < 1.24', lt((arb(1) / 3 + 2 * eps) * a.sqrt() / ln2, arb('1.24')))
c1 = (2 * a + 2 * eps).sqrt() / a.root(4)
chk('c1 = sqrt(2a+2eps) a^{-1/4} < 2.27 (valeur %s)' % c1.str(6), lt(c1, arb('2.27')))
chk('6 c1 / ln2 < 20 (avec la valeur vraie de c1 : %s)' % (6 * c1 / ln2).str(6), lt(6 * c1 / ln2, 20))
chk('log2(4/3) <= 1', lt((arb(4) / 3).log() / ln2, 1))
chk('1.001 + 1.24/sqrt(15) < 4/3', lt(arb('1.001') + arb('1.24') / arb(15).sqrt(), arb(4) / 3))
chk('kappa sqrt(15) < 29 (valeur %s)' % (kappa * arb(15).sqrt()).str(6), lt(kappa * arb(15).sqrt(), 29))
chk('kappa = (pi/ln2) sqrt(8/3)', abs(kappa - pi / ln2 * (arb(8) / 3).sqrt()) < arb('1e-50'))
print('== Corollaire 4.2 / theoreme 4.4')
chk('20 * 2^{1/4} <= 24  (absorption dans e1, non ecrite dans la preuve)', lt(20 * arb(2).root(4), 24))
chk('(2/3 + 1/sqrt2)/kappa < 0.19 (valeur %s)' % ((arb(2) / 3 + 1 / sq2) / kappa).str(6), lt((arb(2) / 3 + 1 / sq2) / kappa, arb('0.19')))
chk('(12 + 1/2)/kappa < 1.7 (valeur %s)' % (arb('12.5') / kappa).str(6), lt(arb('12.5') / kappa, arb('1.7')))
# log2(1+x) <= 4 x^{1/4} pour x >= 0 : f(x)=4x^{1/4}-log2(1+x), f(0)=0, f'(x)>0 car
# x^{-3/4}(1+x) ln2 >= 4*3^{-3/4} ln2 > 1 (minimum en x=3).
chk('log2(1+x) <= 4 x^{1/4} : 4 * 3^{-3/4} * ln2 > 1', lt(1, 4 * arb(3).root(4) ** -3 * ln2))

print('== Forme explicite de la saturation (texte corrige T6)')
chk('(1 + 4 sqrt2/3)/kappa < 0.4 (valeur %s)' % ((1 + 4 * sq2 / 3) / kappa).str(6), lt((1 + 4 * sq2 / 3) / kappa, arb('0.4')))
chk('26 sqrt2/kappa < 5 (valeur %s)' % (26 * sq2 / kappa).str(6), lt(26 * sq2 / kappa, 5))
bcond = (arb(4) / (3 * arb(15).sqrt()) + 26 * arb(4000).root(4) ** -1) / kappa
chk('h>=15m et mh>=4000 => e2/(kappa sqrt(mh)) <= 1/2 (majorant %s)' % bcond.str(6), lt(bcond, arb('0.5')))
# 2(mh)^{-1/2} <= 2(mh)^{-1/4} pour mh >= 1 : absorption analytique, pas de test.

print('\n== Controle numerique (non probant) des quantites modulaires, mpmath 60 chiffres')
mp.mp.dps = 60


def Aq(q):
    t3 = mp.jtheta(3, 0, mp.sqrt(q)); t4 = mp.jtheta(4, 0, mp.sqrt(q))
    P = mp.qp(q)
    return (t3 ** 4 + t4 ** 4) / (2 * P ** 4)


def Bq(q):
    t2 = mp.jtheta(2, 0, mp.sqrt(q))
    return t2 ** 4 / (2 * mp.qp(q) ** 4)
# mpmath : jtheta(3,0,nome) = sum nome^{n^2}; nome = q^{1/2} donne sum q^{n^2/2}.


am = 2 * mp.pi ** 2 / 3
worst_id = 0
for b in [mp.mpf(x) / 20 for x in range(1, 21)] + [mp.mpf(x) for x in (2, 3, 5, 2 * mp.pi, 10, 30)]:
    qt = mp.exp(-4 * mp.pi ** 2 / b)
    lhsA = Aq(mp.exp(-b)); rhsA = mp.exp(am / b - b / 6) / 2 * (Aq(qt) + 3 * Bq(qt))
    lhsB = Bq(mp.exp(-b)); rhsB = mp.exp(am / b - b / 6) / 2 * (Aq(qt) - Bq(qt))
    worst_id = max(worst_id, abs(lhsA / rhsA - 1), abs(lhsB / rhsB - 1))
print('identites modulaires : ecart relatif max = %s' % mp.nstr(worst_id, 3))
chk('identites modulaires (numerique, 60 chiffres)', worst_id < mp.mpf('1e-40'))


def gfun(b):
    y = mp.exp(-2 * mp.pi ** 2 / b)
    return mp.log(Aq(y * y) + 3 * Bq(y * y))


def rfun(b):
    return mp.log(Bq(mp.exp(-b)) / Aq(mp.exp(-b)))


mx = {k: 0 for k in ('g', "g'", "g''", 'r', "r'", "r''")}
sgn_ok = True
for b in [mp.mpf(x) / 40 for x in range(4, 41)]:
    vals = {'g': gfun(b), "g'": mp.diff(gfun, b, 1), "g''": mp.diff(gfun, b, 2),
            'r': rfun(b), "r'": mp.diff(rfun, b, 1), "r''": mp.diff(rfun, b, 2)}
    for kk, vv in vals.items():
        mx[kk] = max(mx[kk], abs(vv))
    tol = mp.mpf(10) ** -50   # bruit de la differentiation numerique
    ok_b = vals['g'] >= -tol and vals["g'"] >= -tol and vals['r'] <= tol and vals["r'"] <= tol
    if not ok_b:
        print('  signe douteux en beta=%s : %s' % (mp.nstr(b, 4), {k: mp.nstr(v, 3) for k, v in vals.items()}))
    sgn_ok &= ok_b
print('max sur beta in [0.1,1] : ' + ', '.join('|%s|=%s' % (k, mp.nstr(v, 3)) for k, v in mx.items()))
chk('signes g,g\' >= 0, r,r\' <= 0 et toutes bornes < 1e-4 (numerique)', sgn_ok and max(mx.values()) < mp.mpf('1e-4'))
worst = mp.inf
for j in range(-170, 301):  # beta in [10^-1.7, 10^3] (qp ne converge plus en dessous)
    b = mp.mpf(10) ** (mp.mpf(j) / 100)
    worst = min(worst, am / b - mp.log(Aq(mp.exp(-b))))
print('min de a/beta - log A(e^-beta) sur beta in [0.02,1000] : %s' % mp.nstr(worst, 5))
chk('A(e^{-beta}) < e^{a/beta} sur la grille (numerique)', worst > 0)

print('\nVERDICT : %d/%d PASS' % (sum(res), len(res)))
sys.exit(0 if all(res) else 1)
