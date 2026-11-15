#!/usr/bin/env python3
"""V2-S4 -- recertification independante (Arb, python-flint) des constantes
decimales utilisees par V2-S3, y compris celles qui different de S2-bis
(b < 0.49 au lieu de b < 1/2). Chaque PASS est une inegalite prouvee par
arithmetique d'intervalles ; un FAIL signale une constante non garantie."""
import sys
from flint import arb, ctx

ctx.prec = 256
R = []


def chk(name, ok):
    R.append(bool(ok))
    print(('PASS ' if ok else 'FAIL ') + name)


def lt(x, y):
    return (arb(y) - arb(x)) > 0


pi, ln2 = arb.pi(), arb(2).log()
a = 2 * pi ** 2 / 3
kappa = 2 * a.sqrt() / ln2
s2 = arb(2).sqrt()
print('kappa =', kappa.str(12))
chk('kappa = (pi/ln2) sqrt(8/3)', abs(kappa - pi / ln2 * (arb(8) / 3).sqrt()) < arb('1e-60'))

# (i)  e1 absorbe 20*(2mh)^{1/4}
chk('20*2^{1/4} < 24', lt(20 * arb(2).root(4), 24))
# (i)  coefficient de sqrt(m/h) : (1/sqrt2)[(4/3)/(sqrt2 kappa) + 1/kappa] = (2/3+1/sqrt2)/kappa
c_i1 = (arb(2) / 3 + 1 / s2) / kappa
chk('(2/3+1/sqrt2)/kappa < 0.19  [%s]' % c_i1.str(6), lt(c_i1, '0.19'))
# (i)  coefficient de (mh)^{-1/4} : 12/kappa + (1/2)/kappa (terme (mh)^{-1/2} absorbe)
c_i2 = arb('12.5') / kappa
chk('12.5/kappa < 1.7  [%s]' % c_i2.str(6), lt(c_i2, '1.7'))

# (iii) b = e2/X, pire cas h=15m, mh=4000 (les deux bornes sont monotones)
b_v2s3 = (arb(4) / 3 / arb(15).sqrt() + 26 / arb(4000).root(4)) / kappa
chk('b <= [(4/3)/sqrt15 + 26*4000^{-1/4}]/kappa < 0.49  [%s]' % b_v2s3.str(6), lt(b_v2s3, '0.49'))
b_exact = (arb(4) / 3 / arb(15).sqrt() + 24 / arb(4000).root(4) + 2 / arb(4000).sqrt()) / kappa
print('     (valeur sans absorption 2/sqrt(mh) <= 2 (mh)^{-1/4} : %s)' % b_exact.str(6))
c_iii1 = (1 + 4 * s2 / 3) / kappa
chk('(1+4sqrt2/3)/kappa < 0.4  [%s]' % c_iii1.str(6), lt(c_iii1, '0.4'))
c_iii2 = 26 * s2 / kappa
chk('26 sqrt2/kappa < 5  [%s]  (marge %s)' % (c_iii2.str(6), (5 - c_iii2).str(3)), lt(c_iii2, 5))
# Corollaire 4.2 S2-bis : log2(1+x) <= 4 x^{1/4}
chk('4*3^{-3/4}*ln2 > 1 (=> log2(1+x) <= 4x^{1/4})', lt(1, 4 * arb(3).root(4) ** -3 * ln2))

# Seuil asymptotique (m -> infini, h = r m) : exposant de Chernoff exact
# phi(r) = -1 + kappa*sqrt(r - 1/6) + o(1); la borne (2) exige 1 + kappa(sqrt r - sqrt(r-1/6)) <= 4/3
for r in (15, 7.5, 3.5):
    rr = arb(r)
    loss = 1 + kappa * (rr.sqrt() - (rr - arb(1) / 6).sqrt())
    print('     perte lineaire par bloc a h/m=%s : %s  (budget 4/3)' % (r, loss.str(6)))
loss15 = 1 + kappa * (arb(15).sqrt() - (arb(15) - arb(1) / 6).sqrt())
chk('perte asymptotique par bloc a h=15m < 4/3  [%s]' % loss15.str(6), lt(loss15, arb(4) / 3))

print('\nVERDICT : %d/%d PASS' % (sum(R), len(R)))
sys.exit(0 if all(R) else 1)
