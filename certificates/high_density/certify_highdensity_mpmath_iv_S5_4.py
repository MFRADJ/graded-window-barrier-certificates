#!/usr/bin/env python3
"""S5.4 independent interval checks of high-density DECIMAL PREMISES.

This is newly written for S5.4 and is NOT the historical S2-bis Arb script.
Each comparison is interval arithmetic using mpmath.iv (80 decimal digits).
It certifies the displayed numerical inequalities *under the explicit
analytic reductions in the manuscript*, not those reductions themselves.
The three historical numerical-only tests of the Arb script are not counted.
"""
from mpmath import iv, __version__ as version
import sys
iv.dps = 80
z = lambda s: iv.mpf(str(s))
res = []

def chk(label, statement):
    is_ok = statement is True
    res.append((label, is_ok))
    print(('PASS ' if is_ok else 'FAIL ') + label)
    if not is_ok:
        raise SystemExit('STOP: interval did not certify ' + label)

pi, ln2 = iv.pi, iv.ln(2)
a = 2 * pi**2 / 3
c = 2*pi**2
kappa = 2*iv.sqrt(a)/ln2
sq2 = iv.sqrt(2)
eps = z('0.0001')
R = z(1)/16
v=R**2
P4 = 1/(1-4*v/(1-v)**2)
Abd = ((1+R)/(1-R))**4*P4
Bbd = 8*R/(1-v)**4*P4
Rs=iv.exp(-pi);qs=iv.exp(-2*pi)
Aqs=(((1+Rs)/(1-Rs))**4 + ((1-Rs)/(1+Rs))**4)/2/(1-4*qs/(1-qs)**2)
Bqs=8*Rs/(1-qs)**4/(1-4*qs/(1-qs)**2)
c1=iv.sqrt(2*a+2*eps)/iv.sqrt(iv.sqrt(a))
bcond=(z(4)/(3*iv.sqrt(15)) +26/iv.sqrt(iv.sqrt(4000)))/kappa

chk('19.7 < 2*pi^2 < 20', z('19.7') < c and c < z(20))
chk('exp(-c) < 3e-9', iv.exp(-c) < z('3e-9'))
chk('c exp(-c) < 6e-8', c*iv.exp(-c) < z('6e-8'))
chk('c^2 exp(-c) < 1.32e-6', c**2 * iv.exp(-c) < z('1.32e-6'))
chk('c > 4', c > 4)
chk('c > 2', c > 2)
chk('bound A on Cauchy circle < 1.68', Abd < z('1.68'))
chk('bound B on Cauchy circle < 0.52', Bbd < z('0.52'))
chk('bound F and G < 4', Abd+3*Bbd < 4)
chk('exp(-pi) < 1/23', Rs < z(1)/23)
chk('exp(-2pi) < 1/500', qs < z(1)/500)
chk('A(exp(-2pi)) bound < 1.08', Aqs < z('1.08'))
chk('B(exp(-2pi)) bound < 0.36', Bqs < z('0.36'))
chk('F at exp(-pi) < 2.2', Aqs+3*Bqs < z('2.2'))
chk('2.2 < 2 exp(1/6)', z('2.2') < 2*iv.exp(z(1)/6))
chk('0.08*2*pi < a', z('0.08')*2*pi < a)
chk('Gibbs mu(1)+2sigma(1) bound < 14.01', a+z(1)/6+eps+2*iv.sqrt(2*a+2*eps) < z('14.01'))
chk('1+2eps/ln2 <1.001', 1+2*eps/ln2 < z('1.001'))
chk('(1/3+2eps)sqrt(a)/ln2 <1.24', (z(1)/3+2*eps)*iv.sqrt(a)/ln2 < z('1.24'))
chk('c1<2.27', c1 < z('2.27'))
chk('6 c1/ln2<20', 6*c1/ln2 < 20)
chk('log2(4/3)<1', iv.ln(z(4)/3)/ln2 < 1)
chk('1.001 +1.24/sqrt15<4/3', z('1.001') + z('1.24')/iv.sqrt(15) < z(4)/3)
chk('kappa sqrt15<29', kappa*iv.sqrt(15) < 29)
chk('20*2^(1/4)<24', 20*iv.sqrt(sq2) < 24)
chk('(2/3 +1/sqrt2)/kappa<0.19', (z(2)/3+1/sq2)/kappa < z('0.19'))
chk('12.5/kappa<1.7', z('12.5')/kappa < z('1.7'))
chk('4*3^(-3/4)*ln2>1', 4*ln2/(iv.sqrt(iv.sqrt(3))**3) > 1)
chk('(1+4sqrt2/3)/kappa<0.4', (1+4*sq2/3)/kappa < z('0.4'))
chk('26sqrt2/kappa<5', 26*sq2/kappa < 5)
chk('b condition<0.5', bcond < z('0.5'))
chk('stronger S4-only b condition<0.49', bcond < z('0.49'))
print('==> V2-S5.4 HIGH-DENSITY DECIMAL PREMISES: %d/%d INTERVAL CHECKS PASSED' % (sum(ok for _,ok in res),len(res)))
print('IMPORTANT: This is new mpmath.iv verification, not a rerun of S2-bis Arb.')
print('Precision: 80 decimal digits; mpmath: ' + version)
